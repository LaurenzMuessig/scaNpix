function eegFilt = lfpFilter(eeg, peakFreq, sampleRate, options)
% lfpFilter - Band pass filter EEG around some peak frequency (e.g. theta)
% Zero phase FIR filter (1s long, Blackman window) around peakFreq +/- filtHalfBandWidth
%
% Usage:
%
%       eegFilt = scanpix.lfp.lfpFilter(eeg, peakFreq, sampleRate)
%       eegFilt = scanpix.lfp.lfpFilter(eeg, peakFreq, sampleRate, 'filtHalfBandWidth', 3)
%
% Inputs:   eeg        - single EEG trace (row or column vector, in uV/V; not int16)
%           peakFreq   - peak frequency around which to filter eeg (in Hz)
%           sampleRate - sample rate of eeg (in Hz)
%
% Options:
%
%   'filtHalfBandWidth'    3,    - filter around peakFreq +/- filtHalfBandWidth
%
% Outputs:  eegFilt    - filtered EEG trace (same orientation as input)
%
% This is based on 'eeg_filter' in Scan (by TW)
% See also: scanpix.lfp.getThetaPhase
%
% TW/LM 2020
%%
arguments
    eeg {mustBeVector, mustBeFloat}
    peakFreq (1,1) {mustBeNumeric}
    sampleRate  (1,1) {mustBeNumeric}
    options.filtHalfBandWidth (1,1) {mustBeNumeric} = 3;
end

%% Filter
passBand       = ([-1 1] .* options.filtHalfBandWidth  +  peakFreq);
window         = fir1( round(sampleRate), passBand./(sampleRate/2), blackman(round(sampleRate)+1) );
eegFilt        = filtfilt(window, 1, eeg );

end
