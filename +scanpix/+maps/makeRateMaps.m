function [ rMaps, sm_pMaps, sm_spkMaps ] = makeRateMaps(obj, trialInd, mapType, addPosFilter, cellInd)
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

%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
    mapType (1,:) {mustBeMember(mapType,{'rate','pos'})} = 'rate';
    addPosFilter {mustBeNumericOrLogical} = false(size(obj.posData.XY{trialInd},1),1);
    cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);
end

%%
% expand smooth kernel
kernel    = ones(obj.mapParams.rate.kernel);
% only relevant for npix data - check for interpolated Fs
if strcmp(obj.type,'npix') && isKey(obj.params,'InterpPos2PosFs') && obj.params('InterpPos2PosFs')
    sampleTimes = [];
else
    sampleTimes = obj.spikeData.sampleT{trialInd};
end

% data from object
positions                 = obj.posData.XY{trialInd};
positions(addPosFilter,:) = NaN;
%
spkTimes                  = obj.spikeData.spk_Times{trialInd}(cellInd);

%% speed filter
if obj.mapParams.rate.speedFilterFlagRMaps 
    speedFilter              = obj.posData.speed{trialInd} <= obj.mapParams.rate.speedFilterLimits(1) | obj.posData.speed{trialInd} > obj.mapParams.rate.speedFilterLimits(2);
    positions(speedFilter,:) = NaN;
end

%% posMap
binSizePix = floor( obj.trialMetaData(trialInd).ppm/100 * obj.mapParams.rate.binSizeSpat ); % this many pix in one bin  %%%%% @TOM: Is 'floor' correct here??
% % bin positions
% % note to comparison to old Scan - to replicate the exact map from the original version you would have to do the following:
% positions = bsxfun(@minus,positions, floor( min(positions,[],1)/binSizePix )*binSizePix ); 
if ~obj.trialMetaData(trialInd).PosIsFitToEnv{1}
    positions = positions - min(positions) + eps; % if you have't fit the positions to the environment and your environment is badly sampled along some edge(s), the resulting ratemap will not be binned correctly
end
% bin
posBinned = fliplr( ceil( positions ./ binSizePix ) ); % swap xy to image coordinates
% get size in bins
if isempty(obj.trialMetaData(trialInd).envSize)
    nBins    = [max(posBinned(:,1)) max(posBinned(:,2))];  % get env size from positions - will be off if env isn't sampled to full extent
else
    envSzPix = obj.trialMetaData(trialInd).envSize ./ 100 .* obj.trialMetaData(trialInd).ppm;
    nBins    = fliplr( ceil( envSzPix ./ binSizePix ) ); %+ min(posBinned)-1; 
end

% raw pos map
posMapRaw    = accumarray(posBinned(~isnan(posBinned(:,1)),:), 1, nBins ) ./ obj.trialMetaData(trialInd).posFs;
unVisPos     = posMapRaw == 0; % keep record of unvisited positions

%% make maps
% pre-allocate
rMaps                     = cell(length(spkTimes),1);
sm_spkMaps                = cell(length(spkTimes),1);
% if using boxcar filtering only need to do pos map once
if strcmp(obj.mapParams.rate.smooth,'boxcar') || strcmp(mapType,'pos')
    visMask               = ones(size(posMapRaw));
    visMask(unVisPos)     = 0;
    visMask_sm            = imfilter(visMask, kernel);
    sm_pMaps{1}           = imfilter(posMapRaw, kernel) ./ visMask_sm;
    sm_pMaps{1}(unVisPos) = NaN;
    %
    if strcmp(mapType,'pos'); return; end
else
    sm_pMaps              = cell(length(spkTimes), 1); % pre-allocate for adaptive smoothing
end


if obj.mapParams.rate.showWaitBar; hWait = waitbar(0); end

for i = 1:length(spkTimes)
    
    if isempty(spkTimes{i})
        rMaps{i}             = zeros(size(posMapRaw));
        rMaps{i}(unVisPos)   = NaN;
        % need a dummy pos map for cells with 0 spikes when using adaptive smooth
        if strcmp(obj.mapParams.rate.smooth,'adaptive')
            tmpPos           = posMapRaw;    
            tmpPos(unVisPos) = NaN;
            sm_pMaps{i}      = tmpPos;
        end
        continue
    end
    
    % spike Map
    if isempty(sampleTimes) 
        spkPosBinInd = ceil(spkTimes{i} .* obj.trialMetaData(trialInd).posFs ); 
    else
        % as sample times in e.g. neuropixel can have some jitter we can't just bin by sample rate
        [~, spkPosBinInd] = arrayfun(@(x) min(abs(sampleTimes - x)), spkTimes{i}, 'UniformOutput', 0); % this is ~2x faster than running min() on whole array at once
        spkPosBinInd      = cell2mat(spkPosBinInd);
        % [~, spkPosBinInd] = min(abs(bsxfun(@minus, sampleTimes, spkTimes{i}.')), [], 1);

    end
    spkPosBinned     = posBinned(spkPosBinInd,:);
    spkMapRaw        = accumarray(spkPosBinned(~isnan(spkPosBinned(:,1)),:), 1, nBins);
    
    % smooth
    switch obj.mapParams.rate.smooth
        case 'boxcar'
            % smooth
            sm_spkMaps{i}                          = imfilter(spkMapRaw,kernel) ./ visMask_sm;
            sm_spkMaps{1}(unVisPos)                = NaN;
            %rate map
            rMaps{i}                               = sm_spkMaps{i} ./ sm_pMaps{1};   
        case 'adaptive'
            [rMaps{i}, sm_spkMaps{i}, sm_pMaps{i}] = scanpix.maps.adaptiveSmooth(posMapRaw,spkMapRaw,obj.mapParams.rate.alpha); % SCAN function  
    end 
    
    if obj.mapParams.rate.showWaitBar; waitbar(i/length(spkTimes),hWait,sprintf('Making those Rate Maps... %i/%i done.',i,length(spkTimes))); end
end
if obj.mapParams.rate.showWaitBar; close(hWait); end

if obj.mapParams.rate.trimNaNs
    rMaps      = trimNaNs(rMaps);
    sm_spkMaps = trimNaNs(sm_spkMaps);
    sm_pMaps   = trimNaNs(sm_pMaps);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maps = trimNaNs(maps)

nanIndC1 = find(sum(isnan(maps{1}),1) ~= size(maps{1},1),1,'first');
nanIndC2 = find(sum(isnan(maps{1}),1) ~= size(maps{1},1),1,'last');
nanIndR1 = find(sum(isnan(maps{1}),2) ~= size(maps{1},2),1,'first');
nanIndR2 = find(sum(isnan(maps{1}),2) ~= size(maps{1},2),1,'last');
maps     = cellfun(@(x) x(nanIndR1:nanIndR2,nanIndC1:nanIndC2),maps,'uni',0);

end


