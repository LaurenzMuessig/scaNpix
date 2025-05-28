function [dirMaps, dirPosMap] = makeDirMaps(obj, trialInd, options)
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

%% 
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
    options.addDirFilter {mustBeNumericOrLogical} = false(size(obj.posData.direction{trialInd},1),1);
    options.cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);
end

%%
% only relevant for npix data - check for interpolated Fs
if strcmp(obj.type,'npix') && isKey(obj.params,'InterpPos2PosFs') && obj.params('InterpPos2PosFs')
    sampleTimes = [];
else
    sampleTimes = obj.spikeData.sampleT{trialInd};
end

% data from object
HDirections                         = obj.posData.direction{trialInd};
HDirections(options.addDirFilter,1) = NaN;
%
spkTimes                            = obj.spikeData.spk_Times{trialInd}(options.cellInd);

%% speed filter
if obj.mapParams.dir.speedFilterFlagDMaps 
    speedFilter                = obj.posData.speed{trialInd} <= obj.mapParams.dir.speedFilterLimitLow | obj.posData.speed{trialInd} > obj.mapParams.dir.speedFilterLimitHigh;
    HDirections(speedFilter,:) = NaN;
end

%% occupancy Map
nBins                   = ceil(360/obj.mapParams.dir.binSizeDir);
HDBinned                = ceil(HDirections ./ obj.mapParams.dir.binSizeDir);
HDBinned(HDBinned == 0) = max(HDBinned); % bin=0 is the same as last bin
occMapRaw               = accumarray(HDBinned(~isnan(HDBinned)),1,[nBins 1]) ./ obj.trialMetaData(trialInd).posFs  ;

%% make maps
% smooth occupancy map
kernel          = ones(obj.mapParams.dir.dirSmoothKern,1) ./ obj.mapParams.dir.dirSmoothKern;
dirPosMap       = imfilter(occMapRaw,kernel,'circular');
 
if obj.mapParams.dir.showWaitBar; hWait = waitbar(0); end

% pre-allocate
dirMaps          = cell(length(spkTimes),1);
for i = 1:length(spkTimes)
    
    if isempty(spkTimes{i})
         dirMaps{i} = zeros(size(dirPosMap));
         continue
    end
    
    % spike Map
    if isempty(sampleTimes)
        spkPosBinInd = ceil(spkTimes{i} .* obj.trialMetaData(trialInd).posFs ); 
    else
        % as sample times can be somewhat irregular we can't just bin by sample rate for e.g. neuropixel data
%         [~, spkPosBinInd] = arrayfun(@(x) min(abs(sampleTimes - x)), spkTimes{i}, 'UniformOutput', 0); % this is ~2x faster than running min() on whole array at once
    [~, spkPosBinInd] = min(abs(bsxfun(@minus, sampleTimes, spkTimes{i}.')), [], 1); % this seems to be most efficient way 
%         spkPosBinInd = cell2mat(spkPosBinInd);
    end
    
    spkPosBinned    = HDBinned(spkPosBinInd,:);
    spkMapRaw       = accumarray(spkPosBinned(~isnan(spkPosBinned(:,1)),:),1,[nBins 1]);
    
    % smooth     
    sm_spkDirMap    = imfilter(spkMapRaw,kernel,'circular');
    % dir map
    dirMaps{i}      = sm_spkDirMap ./ dirPosMap;
        
    if obj.mapParams.dir.showWaitBar; waitbar(i/length(spkTimes),hWait,sprintf('Making those Dir Maps... %i/%i done.',i,length(spkTimes))); end

end

if obj.mapParams.dir.showWaitBar; close(hWait); end

end


