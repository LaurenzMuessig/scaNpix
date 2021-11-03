function [dirMaps, dirPosMap] = makeDirMaps(spkTimes, HDirections, sampleTimes, speed, varargin)
% makeDirMaps - generate directional firing rate maps from spike times 
% and heading direction of animal
%
% package: scanpix.maps
%
% Syntax:
%       scanpix.maps.makeDirMaps(spkTimes, HDirections, sampleTimes)
%       scanpix.maps.makeDirMaps(spkTimes, HDirections, sampleTimes, speed)
%       scanpix.maps.makeDirMaps(spkTimes, HDirections, [], speed)
%       scanpix.maps.makeDirMaps(__, prmsStruct )
%       scanpix.maps.makeDirMaps(__, Name-Value comma separated list )
%
% Inputs:
%    spkTimes    - 1xnCell cell array of spike times in seconds
%    HDirections - numeric array of head directions (degrees an radians are
%                  fine)
%    sampleTimes - sample times in s for each pos sample in neural data
%                  time - can be left empty, e.g. for dacq where we have
%                  excatly 50Hz, but for neuropixel there will be a small
%                  amount of jitter between samples (typically sub ms)
%    speed       - (optional) nPosSamplesx1 array of running speeds in cm/s 
%                  (ommit/leave empty if no speed filtering)              
%    varargin    - prmsStruct: structure with parameter fields to be changed from defaults
%                - name-value: comma separated list of name-value pairs
%
% Outputs:
%   dirMaps      - smoothed directional firing rate map
%
% see also: scanpix.maps.addMaps; scanpix.maps.makeDirMaps; scanpix.maps.makeLinRMaps;
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: CHECK IT'S ALL GOOD

%% params
prms.dirSmoothKern        = 5;         % in bins
prms.binSizeDir           = 6;         % in degrees
prms.speedFilterFlagDMaps = 0;  % y/n
prms.speedFilterLimits    = [2.5 400]; % in cm/s
prms.posFs                = 50;        % in Hz
prms.showWaitBar          = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in case just 1 cell
if ~iscell(spkTimes)
    spkTimes = {spkTimes};
end
% we want degrees
if max(HDirections) <= 2*pi
    HDirections = HDirections * 180/pi;
end

if nargin < 3
    sampleTimes = [];
end

%% speed filter
if nargin > 3 && prms.speedFilterFlagDMaps && ~isempty(speed)
    speedFilter            = speed <= prms.speedFilterLimits(1) | speed > prms.speedFilterLimits(2);
else
    speedFilter            = false(length(HDirections),1); 
end
HDirections(speedFilter,:) = NaN;

%% occupancy Map
nBins                   = ceil(360/prms.binSizeDir);
HDBinned                = ceil(HDirections ./ prms.binSizeDir);
HDBinned(HDBinned == 0) = max(HDBinned); % bin=0 is the same as last bin
occMapRaw               = accumarray(HDBinned(~isnan(HDBinned)),1,[nBins 1]) ./ prms.posFs;

%% make maps
% smooth occupancy map
kernel          = ones(prms.dirSmoothKern,1) ./ prms.dirSmoothKern;
dirPosMap       = imfilter(occMapRaw,kernel,'circular');
 
if prms.showWaitBar; hWait = waitbar(0); end

% pre-allocate
dirMaps          = cell(length(spkTimes),1);
for i = 1:length(spkTimes)
    % spike Map
    if isempty(sampleTimes)
        spkPosBinInd = ceil(spkTimes{i} .* prms.posFs ); 
    else
        % as sample times can be somewhat irregular we can't just bin by sample rate for e.g. neuropixel data
        [~, spkPosBinInd] = arrayfun(@(x) min(abs(sampleTimes - x)), spkTimes{i}, 'UniformOutput', 0); % this is ~2x faster than running min() on whole array at once
        spkPosBinInd = cell2mat(spkPosBinInd);
    end
    
    spkPosBinned    = HDBinned(spkPosBinInd,:);
    spkMapRaw       = accumarray(spkPosBinned(~isnan(spkPosBinned(:,1)),:),1,[nBins 1]);
    
    % smooth     
    sm_spkDirMap    = imfilter(spkMapRaw,kernel,'circular');
    % dir map
    dirMaps{i}      = sm_spkDirMap ./ dirPosMap;
        
    if prms.showWaitBar; waitbar(i/length(spkTimes),hWait,sprintf('Making those Dir Maps... %i/%i done.',i,length(spkTimes))); end

end

if prms.showWaitBar; close(hWait); end

end


