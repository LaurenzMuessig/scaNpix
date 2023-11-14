function [ lin_rMaps, lin_pMaps, linPos ] = makeLinRMaps(spikeTimes,positions,sampleTimes,direction,speed,trackProps,varargin)
% makeLinRMaps - Make linear rate map
% Function will do the position linearisation as well
% package: scanpix.maps
%
%
% Syntax:
%       scanpix.maps.makeLinRMaps(spkTimes, positions, sampleTimes, ppm)
%       scanpix.maps.makeLinRMaps(spkTimes, positions, sampleTimes, ppm, speed)
%       scanpix.maps.makeLinRMaps(spkTimes, positions, [], ppm, speed)
%       scanpix.maps.makeLinRMaps(__, prmsStruct )
%       scanpix.maps.makeLinRMaps(__, Name-Value comma separated list )
%
% Inputs:
%    spkTimes    - 1xnCell cell array of spike times in seconds
%    positions   - numeric array of xy postions
%    sampleTimes - sample times in s for each pos sample in neural data
%                  time - can be left empty, e.g. for dacq where we have
%                  excatly 50Hz, but for neuropixel there will be a small
%                  amount of jitter between samples (typically sub ms)
%    direction   - numeric array of head directions in degrees
%    speed       - (optional) nPosSamplesx1 array of running speeds in cm/s
%                  (ommit/leave empty if no speed filtering)  
%    trackProps  - struct with track propereties. needs to include fields 'type' ('sqtrack' or 'lintrack'), 
%                  'length' (length in cm; for sq track length of 1 arm), 'ppm' (pix/m from tracking) 
%                  and 'posFs' (position sample rate)
%    varargin    - prmsStruct: structure with parameter fields to be changed from defaults
%                - name-value: comma separated list of name-value pairs
%
% Outputs:
%   lin_rMaps    - smoothed linear firing rate maps {full, CW, CCW} 
%   lin_pMaps    - smoothed linear position maps {full, CW, CCW}
%   linPos       - array of linear position
%   lin_rMaps_normed - normalised and sorted by peak position on track
%                      version of lin_rMaps - {full, sortIndex - full; CW,
%                      sortIndex - CW; CCW, sortIndex - CCW}
%
% see also: scanpix.maps.addMaps; scanpix.maps.linearisePosData;
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% params
% Parameters for linearisation of position data 
prms.minDwellForEdge      = 1; % in s
prms.durThrCohRun         = 2; % In seconds
prms.filtSigmaForRunDir   = 3;  % Units=sec. Sigma of the Gaussian filter to pre-filter the data before finding CW and CCW runs. Kernel is 2*sigma in length.
prms.durThrJump           = 2;  % In seconds
prms.gradThrForJumpSmooth = 2.5;
prms.runDimension         = [];                   % This is the dimension to run (X=1, Y=2)
prms.dirTolerance         = 70;  % Tolerance for heading direction, relative to arm direction, for calculating run direction (degrees)

% linear map params
prms.showWaitBar          = false;
prms.binSizeLinMaps       = 2.5; % bin size (cm^2)
prms.smoothFlagLinMaps    = 1;  % y/n
prms.smoothKernelSD       = 2;  % SD, in bins
prms.speedFilterFlagLMaps = 1;  % y/n
prms.speedFilterLimits    = [2.5 400]; % See setting above - lock settings for linear and open field together.
prms.posIsCircular        = 0;
prms.remTrackEnds         = 0;  % Remove this many bins from each end of the linear track. Set to 0 to remove none.
prms.normSort             = 1;

% ---------------------------------------------------------------------------------- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %

% convert to radians
prms.dirTolerance = prms.dirTolerance * pi/180; 

if ~isempty(regexp(lower(trackProps.type),'sq','once'))
    prms.posIsCircular = 1;
end

%% linearise position
[linPos, dirInd ] = scanpix.maps.linearisePosData(positions,direction,trackProps, prms);
dirInd = [dirInd == 1, dirInd == 2]; % make logical

linPos = repmat(linPos,1,3);
linPos(~dirInd(:,1),2) = NaN;
linPos(~dirInd(:,2),3) = NaN;

%% make rate maps

% Make smoothing kernel %
kernel = fspecial('Gaussian',[1 5*prms.smoothKernelSD], prms.smoothKernelSD); % should kernel be shorter??

% speed filter
if prms.speedFilterFlagLMaps && ~isempty(speed)
    speedInd = speed <= prms.speedFilterLimits(1) | speed > prms.speedFilterLimits(2);
    linPos(speedInd,:) = NaN;
end
  
% bin pos
binSizePix   = floor( (trackProps.ppm/100) * prms.binSizeLinMaps ); % is floor the right thing here? 
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
lin_pMaps(1,:) = histcounts(linPosBinned,binList) ./ trackProps.posFs ;
lin_pMaps(2,:) = histcounts(linPosBinned(dirInd(:,1)),binList) ./ trackProps.posFs ;
lin_pMaps(3,:) = histcounts(linPosBinned(dirInd(:,2)),binList) ./ trackProps.posFs ;

% Smooth pos maps %
if prms.smoothFlagLinMaps    
    if strcmpi(trackProps.type,'sqtrack')
        lin_pMaps     = imfilter(lin_pMaps, kernel, 'conv', 'circular');
    else
        smMap         = imfilter(lin_pMaps, kernel, 'conv', 0);
        smKernel      = imfilter(ones(size(lin_pMaps)), kernel, 'conv', 0);
        lin_pMaps     = smMap./smKernel;
    end 
end

if prms.showWaitBar; hWait = waitbar(0); end

% Spike (and therefore rate) maps %
% pre-allocate
lin_rMaps   = cell(3,1);
% nSpkPostExc = nan( length(spikeTimes), 1 );
for s = 1:length(spikeTimes)
    
    % Spike maps (spike counts) %
    if ~isempty(spikeTimes{s})
        spkMaps = nan(size(lin_pMaps));
        % Convert spike times to pos bin values %
%         if isempty(sampleTimes)
%             temp_ST           = ceil(spikeTimes{s} .* trackProps.posFs ); 
%         else
%             temp_ST = arrayfun(@(x) find(sampleTimes - x > 0,1,'first'), spikeTimes{s}, 'UniformOutput', 0); % as sample times can be somewhat irregular we can't just bin by sample rate
%             temp_ST = cell2mat(temp_ST);
%         end
        if isempty(sampleTimes)
            temp_ST           = ceil(spikeTimes{s} .* trackProps.posFs );
        else
            % as sample times in e.g. neuropixel can have some jitter we can't just bin by sample rate
%             [~, temp_ST] = arrayfun(@(x) min(abs(sampleTimes - x)), spikeTimes{s}, 'UniformOutput', 0); % this is ~2x faster than running min() on whole array at once
%             temp_ST = cell2mat(temp_ST);
            [~, temp_ST] = min(abs(bsxfun(@minus, sampleTimes, spikeTimes{s}.')), [], 1);
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
    
    if prms.smoothFlagLinMaps
        if strcmpi(trackProps.type,'sqtrack')
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
    
    if prms.showWaitBar; waitbar(i/length(spikeTimes),hWait,sprintf('Making those Rate Maps... %i/%i done.',s,length(spkTimes))); end
end

if prms.remTrackEnds > 0
    mapBins = 1:length(lin_rMaps{1});
    set2nanInd = mapBins <= prms.remTrackEnds | mapBins > mapBins(end) - prms.remTrackEnds;
    for i = 1:3
        lin_rMaps{i}(:,set2nanInd) = NaN;
    end
    lin_pMaps(:, set2nanInd) = NaN;
end

if prms.showWaitBar; close(hWait); end

end

