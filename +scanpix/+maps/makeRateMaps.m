function [ rMaps, sm_pMaps, sm_spkMaps ] = makeRateMaps(spkTimes, positions, sampleTimes, ppm, speed, varargin)
% makeRateMaps - Make a spatial rate map
% package: scanpix.maps
%
%
% Syntax:
%       scanpix.maps.makeRateMaps(spkTimes, positions, sampleTimes, ppm)
%       scanpix.maps.makeRateMaps(spkTimes, positions, sampleTimes, ppm, speed)
%       scanpix.maps.makeRateMaps(spkTimes, positions, [], ppm, speed)
%       scanpix.maps.makeRateMaps(__, prmsStruct )
%       scanpix.maps.makeRateMaps(__, Name-Value comma separated list )
%
% Inputs:
%    spkTimes    - 1xnCell cell array of spike times in seconds
%    positions   - numeric array of xy postions
%    sampleTimes - sample times in s for each pos sample in neural data
%                  time - can be left empty, e.g. for dacq where we have
%                  excatly 50Hz, but for neuropixel there will be a small
%                  amount of jitter between samples (typically sub ms)
%    ppm         - pixel/m value from tracking
%    speed       - (optional) nPosSamplesx1 array of running speeds in cm/s
%                  (ommit/leave empty if no speed filtering)  
%    varargin    - prmsStruct: structure with parameter fields to be changed from defaults
%                - name-value: comma separated list of name-value pairs
%
% Outputs:
%   rMaps        - smoothed spatial firing rate map
%   sm_pMaps     - smoothed spatial position map
%   sm_spkMaps   - smoothed spatial spike rate map
%
% see also: scanpix.maps.addMaps; scanpix.maps.makeDirMaps; scanpix.maps.makeLinRMaps;
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: mapsize defined by input rather than max bin visited

%% params
prms.smooth               = 'boxcar'; % 'boxcar; 'adaptive'
prms.smoothKernel         = 5;
prms.binSizeSpat          = 2.5; % in cm^2
prms.PosFs                = 50;
prms.alpha                = 200;
prms.speedFilterFlagRMaps = 0;  % y/n
prms.speedFilterLimits    = [2.5 400];
prms.mapSize              = [];


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
% expand smooth kernel
kernel = ones(prms.smoothKernel);

%% speed filter
if nargin > 4 && prms.speedFilterFlagRMaps && ~isempty(speed)
    speedFilter = speed <= prms.speedFilterLimits(1) | speed > prms.speedFilterLimits(2);
else
    speedFilter = false(length(positions),1); 
end
positions(speedFilter,:) = NaN;

%% posMap
binSizePix = floor( ppm/100 * prms.binSizeSpat ); % this many pix in one bin  %%%%% @TOM: Is 'floor' correct here??
% bin positions
% first we want to shift positions so that the lower edges will fall into bin 1 in both dimensions (in dacq coords that would be dependent on
% window size). So we find the nearest multiple of 'binSizePix' to lower bound of coordinates and subtract that to avoid moving pixels across bin
% edges - this is mainly to keep format backward compatible to the old Scan (not sure this is really necessary)
positions = bsxfun(@minus,positions, floor( min(positions,[],1)/binSizePix )*binSizePix ); 
% if any(min(positions,[],1) == 0)
%     positions(:,min(positions,[],1) == 0) = positions(:,min(positions,[],1) == 0) + 1;
% end
positions(positions == 0) = 1;

posBinned = fliplr( ceil( positions ./ binSizePix ) ); % swap xy to image coordinates (not necesarry, but keeps old Scan format)

if isempty(prms.mapSize)
    nBins = [nanmax(posBinned(:,1)) nanmax(posBinned(:,2))]; % 'guess' map size from data (fine if sampled well)
else
    nBins = prms.mapSize;
end
% raw pos map
posMapRaw  = accumarray(posBinned(~isnan(posBinned(:,1)),:), 1, nBins ) ./ prms.PosFs;
unVisPos   = posMapRaw == 0; % keep record of unvisited positions

%% make maps
% pre-allocate
rMaps                     = cell(length(spkTimes),1);
sm_spkMaps                = cell(length(spkTimes),1);
% if using boxcar filtering only need to do pos map once
if strcmp(prms.smooth,'boxcar')
    visMask               = ones(size(posMapRaw));
    visMask(unVisPos)     = 0;
    visMask_sm            = imfilter(visMask, kernel);
    sm_pMaps{1}           = imfilter(posMapRaw, kernel) ./ visMask_sm;
    sm_pMaps{1}(unVisPos) = NaN;
else
    sm_pMaps              = cell(length(spkTimes), 1); % pre-allocate for adaptive smoothing
end

for i = 1:length(spkTimes)
    % spike Map
    if isempty(sampleTimes)
        spkPosBinInd = ceil(spkTimes{i} .* prms.PosFs ); 
    else
        % as sample times in e.g. neuropixel can have some jitter we can't just bin by sample rate
        spkPosBinInd = arrayfun(@(x) find(sampleTimes - x > 0,1,'first'), spkTimes{i}, 'UniformOutput', 0); 
        spkPosBinInd = cell2mat(spkPosBinInd);
    end
    try
    spkPosBinned     = posBinned(spkPosBinInd,:);
    catch
        t=1;
    end
    spkMapRaw        = accumarray(spkPosBinned(~isnan(spkPosBinned(:,1)),:), 1, nBins);
    
    % smooth
    switch prms.smooth
        case 'boxcar'
            % smooth
            sm_spkMaps{i}                          = imfilter(spkMapRaw,kernel) ./ visMask_sm;
            sm_spkMaps{1}(unVisPos)                = NaN;
            %rate map
            rMaps{i}                               = sm_spkMaps{i} ./ sm_pMaps{1};   
        case 'adaptive'
            [rMaps{i}, sm_spkMaps{i}, sm_pMaps{i}] = scanpix.maps.adaptiveSmooth(posMapRaw,spkMapRaw,prms.alpha); % SCAN function  
        otherwise
            error([prms.smooth ' is not a valid option for smoothing rate maps']);
    end 
end

end

