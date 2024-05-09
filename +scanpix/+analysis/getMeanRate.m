function  meanRates = getMeanRate(spikeTimes,sampleRate,speed,speedLimits)
% getMeanRate - get true mean firing rate while accounting for possible speed filtering 
% package: scanpix.analysis
%
%
% Syntax:
%       meanRates = scanpix.analysis.getMeanRate(spikeTimes,sampleRate,speed)
%       meanRates = scanpix.analysis.getMeanRate(spikeTimes,sampleRate,speed,speedLimits)
%
% Inputs:
%    spikeTimes  - cell array of spike times for a given recording trial
%    sampleRate  - position sample rate
%    speed       - vector of running speed for a given trial (in pos samples)
%    speedLimits - (optional); either 'none' or limits for speed filtering [low high]; if ommitted will default to [2.5 400]
%
%
% Outputs: 
%
%    meanRates    - mean firing rates
%
% LM 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input checks
if nargin < 4
    speedLimits = [2.5 400]; % defaults
end
%
if strcmp(speedLimits,'none')
    speedLimits = [0 10000]; % set arbitrarily high
end
%
if ~iscell(spikeTimes)
    spikeTimes = {spikeTimes};
end

%% get mean rate per cell
speedFilter = speed > speedLimits(1) & speed < speedLimits(2);

if all(speedFilter)
    % as per normal: nSpikes / duration
    meanRates = cell2mat(cellfun(@(x) length(x),spikeTimes,'uni',0)) ./ (length(speedFilter) / sampleRate);
else
    % for speed filtered maps, we need to account for spikes/samples that get removed
    meanRates = cell2mat(cellfun(@(x) sum(speedFilter(ceil(x.*sampleRate)),'omitnan'),spikeTimes,'uni',0)) ./ (sum(speedFilter) / sampleRate);
end

end