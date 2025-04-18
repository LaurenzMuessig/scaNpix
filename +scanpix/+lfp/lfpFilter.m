function [eegFilt, eegPhase, cycleN] = lfpFilter(eeg, peakFreq, sampleRate, options)
% eegFilter - Filter EEG within some frequency band and also extract phase
% We'll try and remove dodgy bits from data by a) making sure that
% individual cycles have > threshold power and b) have lengths within limits
% of 'peakFreq+/-prms.filtHalfBandWidth'.
%
% Usage:
%
%       [eegFilt, eegPhase, cycleN] = eegFilter(eeg, peakFreq)
%       [eegFilt, eegPhase, cycleN] = eegFilter(eeg, peakFreq, optionalInputStruct );
%       [eegFilt, eegPhase, cycleN] = eegFilter(eeg, peakFreq, 'inputName', inputVal, .. etc .. );
%
% Inputs:   eeg      - array of single EEG trace 
%           peakFreq - peak frequency around which to filter eeg (in Hz)
%           varargin - prms.Struct or Name-value pair list
%
% Options:
%
%   'Fs'                   250,  - sample rate for EEG
%   'filtHalfBandWidth'    3,    - filter around peakFreq +/- prms.filtHalfBandWidth
%   'powerThresh'          5,    - threshold (percentile) for min. power/cycle 
%
% Outputs:  eegFilt    - filtered EEG trace
%           eegPhase   - phase in radians of filtered EEG
%           cycleN     - numeric index for cycle number
%
% This is based on 'eeg_filter' in Scan (by TW), but will also get phase and 
% find indices for all cycles 
%
% TW/LM 2020
%%
arguments
    eeg {mustBeNumeric}
    peakFreq (1,1) {mustBeNumeric}
    sampleRate  (1,1) {mustBeNumeric}
    options.filtHalfBandWidth (1,1) {mustBeNumeric} = 3;
    options.powerThresh (1,1) {mustBeNumeric} = 5;
end

%% Filter
% Filter EEG %
passBand       = ([-1 1] .* options.filtHalfBandWidth  +  peakFreq);
window         = fir1( round(sampleRate), passBand./(sampleRate/2), blackman(round(sampleRate)+1) );
eegFilt        = filtfilt(window, 1, eeg );

%% Phase & Cycles
% get phase as well if desired
if nargout > 1
    eegPhase = angle( hilbert(eegFilt) ); % hilbert transform
    eegPhase = mod(eegPhase, 2*pi);       % By wrapping into the range 0 - 2pi, we get the 'classic' theta convention the that oscillation starts at the peak.
    
    % phase transitions, i.e. cycle boundaries
    phaseTrans = diff(eegPhase) < -pi;
    phaseTrans = [true, phaseTrans(1:end-1), true];
    % remove phase slips
    phaseSlips             = [false, diff(unwrap(eegPhase(2:end))) < 0]; % first sample can't be phaseslip
    phaseSlips(end+1)      = phaseSlips(end);
    phaseTrans(phaseSlips) = false;
    % numeric cycle index
    cycleN                 = cumsum(phaseTrans);
    
    % eeg power
    power          = eegFilt.^2;
    powerPerCycle  = accumarray(cycleN(~phaseSlips)', power(~phaseSlips)') ./ accumarray(cycleN(~phaseSlips)',1); % ignore phase slips
    
    % remove bad data
    % power threshold - NEED TO PLAY AROUND WITH THRESHOLD BUT 5th PRCTLE SEEMS LIKE A REASONABLE CHOICE
    thresh         = prctile(powerPerCycle, options.powerThresh);
    badPowerCycle  = find(powerPerCycle < thresh);
    badPowerInd    = ismember(cycleN, badPowerCycle);
    % length threshold
    cycleLength    = accumarray(cycleN',1); % need to include phase slip samples here to get true length of cycles
    cycleLengthLim = ceil(1 ./ passBand * sampleRate); 
    badLengthCycle = find(cycleLength < cycleLengthLim(2) | cycleLength > cycleLengthLim(1));
    badLengthInd   = ismember(cycleN, badLengthCycle);
    
    % remove from data
    cycleN(badLengthInd | badPowerInd)   = NaN;
    eegPhase(badLengthInd | badPowerInd) = NaN;
    
end



end

