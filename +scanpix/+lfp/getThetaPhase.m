function [eegPhase, cycleN, cycleInfo] = getThetaPhase(eegFilt, peakFreq, sampleRate, options)
% getThetaPhase - get phase for each sample of a (theta) filtered EEG and split it into
% individual cycles. We'll try and remove dodgy bits from data by
% a) making sure that individual cycles have enough power (relative to median power/cycle)
% b) cycles have lengths within limits of 'peakFreq +/- filtHalfBandWidth'
% c) removing phase slips (samples where phase runs backwards)
% d) removing partial cycles at start/end of data and at edges of speed filtered periods
%
% Usage:
%
%       [eegPhase, cycleN, cycleInfo] = scanpix.lfp.getThetaPhase(eegFilt, peakFreq, sampleRate)
%       [eegPhase, cycleN, cycleInfo] = scanpix.lfp.getThetaPhase(eegFilt, peakFreq, sampleRate, 'speedFilter', speedFilt)
%       [eegPhase, cycleN, cycleInfo] = scanpix.lfp.getThetaPhase(eegFilt, peakFreq, sampleRate, 'inputName', inputVal, .. etc .. )
%
% Inputs:   eegFilt    - filtered EEG trace (row or column vector), e.g. from scanpix.lfp.lfpFilter
%           peakFreq   - peak frequency eeg was filtered around (in Hz)
%           sampleRate - sample rate of eeg (in Hz)
%
% Options:
%
%   'speedFilter'          [],   - logical index into eeg (true = valid sample), e.g. from
%                                  scanpix.lfp.speedFilter2lfp. Cycles with fewer valid samples
%                                  than 'minRunFrac' are removed and the median power is only
%                                  taken from the remaining cycles. Invalid samples always get
%                                  phase = NaN. Empty = all samples valid
%   'minRunFrac'           1,    - min. fraction of samples per cycle that need to be valid as per
%                                  'speedFilter' (1 = whole cycle has to be valid; e.g. 0.5 =
%                                  keep cycles that are mostly valid)
%   'filtHalfBandWidth'    3,    - half band width eeg was filtered with (sets limits for cycle length)
%   'powerThresh'          0.25, - min. power/cycle as fraction of median power/cycle
%
% Outputs:  eegPhase   - phase in radians of filtered EEG [0 2pi), 0 = peak; NaN for removed data
%           cycleN     - numeric index for cycle number; NaN for removed cycles
%           cycleInfo  - struct with per cycle info (power, length, good/bad) and median power
%
% Phase convention: oscillation starts at the peak, i.e. peak = 0, trough = pi
%
% See also: scanpix.lfp.lfpFilter, scanpix.lfp.speedFilter2lfp
%
% TW/LM 2020 (split off from lfpFilter 2026)
%%
arguments
    eegFilt {mustBeVector, mustBeFloat}
    peakFreq (1,1) {mustBePositive}
    sampleRate (1,1) {mustBePositive}
    options.speedFilter {mustBeNumericOrLogical} = []
    options.minRunFrac (1,1) {mustBeGreaterThan(options.minRunFrac,0), mustBeLessThanOrEqual(options.minRunFrac,1)} = 1
    options.filtHalfBandWidth (1,1) {mustBePositive} = 3
    options.powerThresh (1,1) {mustBeNonnegative} = 0.25
end

if options.filtHalfBandWidth >= peakFreq
    error('scaNpix::lfp::getThetaPhase:''filtHalfBandWidth'' needs to be smaller than ''peakFreq''.');
end

isRowIn = isrow(eegFilt);
eegFilt = eegFilt(:);
nSamp   = length(eegFilt);

if isempty(options.speedFilter)
    speedFilter = true(nSamp,1);
else
    if numel(options.speedFilter) ~= nSamp
        error('scaNpix::lfp::getThetaPhase:''speedFilter'' needs to have same number of samples as eeg (%d vs %d).', numel(options.speedFilter), nSamp);
    end
    speedFilter = logical(options.speedFilter(:));
end

%% phase
eegPhase = angle( hilbert(eegFilt) ); % hilbert transform
eegPhase = mod(eegPhase, 2*pi);       % By wrapping into the range 0 - 2pi, we get the 'classic' theta convention that the oscillation starts at the peak.

%% cycles
% use running max of unwrapped phase: new cycle starts when phase passes 2pi*k (i.e. the peak) for the 1st time, so
% jitter/phase slips around the peak can't create extra cycle boundaries
unwrPhase  = unwrap(eegPhase);
maxPhase   = cummax(unwrPhase);
phaseSlips = unwrPhase < maxPhase; % phase is running backwards (or hasn't caught up again after doing so)
cycleN     = floor(maxPhase ./ (2*pi));
cycleN     = cycleN - cycleN(1) + 1;
nCycles    = cycleN(end);

%% cycle properties
cycleLength   = accumarray(cycleN, 1, [nCycles 1]);                               % includes phase slip samples to get true length of cycles
cycleRunFrac  = accumarray(cycleN, speedFilter, [nCycles 1]) ./ cycleLength;       % fraction of valid samples per cycle
cycleSpeedOK  = cycleRunFrac >= options.minRunFrac;
powerPerCycle = accumarray(cycleN(~phaseSlips), eegFilt(~phaseSlips).^2, [nCycles 1]) ./ accumarray(cycleN(~phaseSlips), 1, [nCycles 1]); % ignore phase slips; NaN if no valid samples

% length limits
passBand       = [-1 1] .* options.filtHalfBandWidth + peakFreq;
cycleLengthOK  = cycleLength >= sampleRate / passBand(2) & cycleLength <= sampleRate / passBand(1);
% partial cycles at start and end of data
edgeCycle      = false(nCycles,1);
edgeCycle([1 nCycles]) = true;

% power threshold - relative to median power of all otherwise valid cycles
candidateCycle = cycleLengthOK & cycleSpeedOK & ~edgeCycle & ~isnan(powerPerCycle);
medianPower    = median(powerPerCycle(candidateCycle));
goodCycle      = candidateCycle & powerPerCycle >= options.powerThresh * medianPower;

%% remove bad data
badCycleSamp           = ~goodCycle(cycleN);
cycleN                 = double(cycleN);
cycleN(badCycleSamp)   = NaN;
eegPhase(badCycleSamp | phaseSlips | ~speedFilter) = NaN; % slip samples keep their cycle ID, but phase is unreliable; same for invalid samples in kept cycles (minRunFrac < 1)

%% output
if isRowIn
    eegPhase = eegPhase';
    cycleN   = cycleN';
end

cycleInfo.power       = powerPerCycle;
cycleInfo.length      = cycleLength;
cycleInfo.lengthOK    = cycleLengthOK;
cycleInfo.runFrac     = cycleRunFrac;
cycleInfo.speedOK     = cycleSpeedOK;
cycleInfo.good        = goodCycle;
cycleInfo.medianPower = medianPower;

end
