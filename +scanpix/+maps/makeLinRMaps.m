function [ lin_rMaps, lin_pMaps ] = makeLinRMaps(obj, trialInd, options )
% makeLinRMaps - Make linear rate map
% Function will do the position linearisation as well
% package: scanpix.maps
%
%
% Syntax:
%       scanpix.maps.makeLinRMaps(spkTimes, positions, sampleTimes, ppm)
%       scanpix.maps.makeLinRMaps(spkTimes, positions, sampleTimes, ppm, speed)
%
% Inputs:
%    spkTimes    - 
%
% Outputs:
%   lin_rMaps    - smoothed linear firing rate maps {full, CW, CCW} 
%   lin_pMaps    - smoothed linear position maps {full, CW, CCW}
%
% see also: scanpix.maps.addMaps; scanpix.maps.linearisePosData;
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% params
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
    options.addPosFilter {mustBeNumericOrLogical} = false(size(obj.posData.XY{trialInd},1),1);
    options.cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);
end

%%
% only relevant for npix data - check for interpolated Fs
if strcmp(obj.type,'npix') && isKey(obj.params,'InterpPos2PosFs') && obj.params('InterpPos2PosFs')
    sampleTimes = [];
else
    sampleTimes = obj.spikeData.sampleT{trialInd};
end

%%
% data from object
linPos                         = obj.posData.linXY{trialInd};
%
if isempty(linPos); warning('scaNpix::maps::makeLinRMaps: You need to linearise the pos data before you ask for lineraised rate maps. Duh!'); return; end
%
linPos(options.addPosFilter,:) = NaN;
%
spikeTimes                     = obj.spikeData.spk_Times{trialInd}(options.cellInd);

%% make rate maps
% Make smoothing kernel %
kernel = fspecial('Gaussian',[1 5*obj.mapParams.lin.smoothKernelSD], obj.mapParams.lin.smoothKernelSD); % should kernel be shorter??

% speed filter
if obj.mapParams.lin.speedFilterFlagLMaps
    speedInd           = obj.posData.speed{trialInd} <= obj.mapParams.lin.speedFilterLimitLow | obj.posData.speed{trialInd} > obj.mapParams.lin.speedFilterLimitHigh;
    linPos(speedInd,:) = NaN;
end
  
% bin pos
binSizePix   = floor( (obj.trialMetaData(trialInd).ppm/100) * obj.mapParams.lin.binSizeLinMaps ); % is floor the right thing here? 
linPosBinned = ceil(linPos(:,1)./binSizePix);
nBins        = max(linPosBinned);  % 
binList      = 0 : nBins;

% I don't think this is a wise thing to do as this means we would need to
% carry over how much was removed for any further analysis - see addition
% in lines 160-67
% If requested, remove ends from all maps (this is to remove the ends of the linear track) %
% Rather than simply make the maps and then chop the ends off, we mark the end positions as 
% nan, then subtract the end zone length to 'shorten' the resulting map. This is so that we
% can (later in the code) count how many spikes are in the end and non-end zones, and thus
% we can remove entirely from the analysis cells with insufficient cells in the non-end zone.
% if prms.remTrackEnds > 0
%     remInd               = linPosBinned <= prms.remTrackEnds | linPosBinned > (nBins - prms.remTrackEnds);
%     linPosBinned(remInd) = nan;
%     linPosBinned         = linPosBinned - prms.remTrackEnds;
%     nBins                = nBins - (prms.remTrackEnds * 2);
% end

% Make position maps %
lin_pMaps      = nan(3, length(binList) - 1);
dirInd         = ~isnan(linPos(:,2:end));      
%
lin_pMaps(1,:) = histcounts(linPosBinned,binList) ./ obj.trialMetaData(trialInd).posFs ;
lin_pMaps(2,:) = histcounts(linPosBinned(dirInd(:,1)),binList) ./ obj.trialMetaData(trialInd).posFs;
lin_pMaps(3,:) = histcounts(linPosBinned(dirInd(:,2)),binList) ./ obj.trialMetaData(trialInd).posFs;

% Smooth pos maps %
if obj.mapParams.lin.smoothFlagLinMaps    
    if ~isempty(regexp(lower(obj.trialMetaData(trialInd).trialType),'sq','once'))
        lin_pMaps     = imfilter(lin_pMaps, kernel, 'conv', 'circular');
    else
        smMap         = imfilter(lin_pMaps, kernel, 'conv', 0);
        smKernel      = imfilter(ones(size(lin_pMaps)), kernel, 'conv', 0);
        lin_pMaps     = smMap./smKernel;
    end 
end

if obj.mapParams.lin.showWaitBar; hWait = waitbar(0); end

% Spike (and therefore rate) maps %
% pre-allocate
lin_rMaps   = cell(3,1);
% nSpkPostExc = nan( length(spikeTimes), 1 );
for s = 1:length(spikeTimes)
    
    % Spike maps (spike counts) %
    if ~isempty(spikeTimes{s})
        spkMaps = nan(size(lin_pMaps));
        %
        if isempty(sampleTimes)
            temp_ST           = ceil(spikeTimes{s} .* obj.trialMetaData(trialInd).posFs );
        else
            % as sample times in e.g. neuropixel can have some jitter we can't just bin by sample rateoptions.
            [~, temp_ST] = arrayfun(@(x) min(abs(sampleTimes - x)), spikeTimes{s}, 'UniformOutput', 0); % this is ~2x faster than running min() on whole array at once
            temp_ST = cell2mat(temp_ST);
            % [~, temp_ST] = min(abs(bsxfun(@minus, sampleTimes, spikeTimes{s}.')), [], 1);
        end
        
        temp_ST(temp_ST == 0) = 1; % this shouldn't happen. not sure why this is here -- CHECK
        % Bin both directions %
        spkLocAll             = linPosBinned(temp_ST,1); %spike locations
        spkMaps(1,:)          = histcounts(spkLocAll,binList);
        % Count spikes from excluded pos samples (so as to get 'true' spike count after exclusions) %
        % This applies to spikes for both directions, is not specific to direction filters.         %
%         nSpkPostExc(s)        = sum( ~isnan( spkLocAll ) );
        % Bin CW runs %
        indST_CW              = temp_ST( dirInd(temp_ST, 1) );
        spkLocCW              = linPosBinned(indST_CW); %spike locations
        spkMaps(2,:)          = histcounts(spkLocCW,binList);
        % Bin CCW runs %
        indST_CCW             = temp_ST(dirInd(temp_ST, 2));
        spkLocCCW             = linPosBinned(indST_CCW); %spike locations
        spkMaps(3,:)          = histcounts(spkLocCCW,binList);
    else
        spkMaps               = zeros(3,nBins);
%         nSpkPostExc(s)        = 0;
    end
    
    if obj.mapParams.lin.smoothFlagLinMaps
        if ~isempty(regexp(lower(obj.trialMetaData(trialInd).trialType),'sq','once'))
            spkMaps = imfilter(spkMaps, kernel, 'conv', 'circular');
        else
            smMap   = imfilter(spkMaps, kernel, 'conv', 0);
            spkMaps = smMap ./ smKernel;
        end
    end
    
    % Rate maps %
    for i = 1:3
        lin_rMaps{i}(s,:) = spkMaps(i,:) ./ lin_pMaps(i,:); %rate
    end
    
    if obj.mapParams.lin.showWaitBar; waitbar(s/length(spikeTimes),hWait,sprintf('Making those Rate Maps... %i/%i done.',s,length(spikeTimes))); end
end

if obj.mapParams.lin.remTrackEnds > 0
    mapBins = 1:length(lin_rMaps{1});
    set2nanInd = mapBins <= obj.mapParams.lin.remTrackEnds | mapBins > mapBins(end) - obj.mapParams.lin.remTrackEnds;
    for i = 1:3
        lin_rMaps{i}(:,set2nanInd) = NaN;
    end
    lin_pMaps(:, set2nanInd) = NaN;
end

if obj.mapParams.lin.showWaitBar; close(hWait); end

end

