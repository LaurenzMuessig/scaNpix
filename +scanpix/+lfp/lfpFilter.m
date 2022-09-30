function [eegFilt, eegPhase, cycleN] = lfpFilter(eeg, peakFreq, varargin)
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

%% Params
prms.Fs                = 250; % sample rate for eeg in Hz
prms.filtHalfBandWidth = 3;   % filter around peakFreq +/- prms.filtHalfBandWidth
prms.powerThresh       = 5;   % percentile


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

%% Filter
% Filter EEG %
passBand       = ([-1 1] .* prms.filtHalfBandWidth  +  peakFreq);
window         = fir1( round(prms.Fs), passBand./(prms.Fs/2), blackman(round(prms.Fs)+1) );
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
    thresh         = prctile(powerPerCycle, prms.powerThresh);
    badPowerCycle  = find(powerPerCycle < thresh);
    badPowerInd    = ismember(cycleN, badPowerCycle);
    % length threshold
    cycleLength    = accumarray(cycleN',1); % need to include phase slip samples here to get true length of cycles
    cycleLengthLim = ceil(1 ./ passBand * prms.Fs); 
    badLengthCycle = find(cycleLength < cycleLengthLim(2) | cycleLength > cycleLengthLim(1));
    badLengthInd   = ismember(cycleN, badLengthCycle);
    
    % remove from data
    cycleN(badLengthInd | badPowerInd)   = NaN;
    eegPhase(badLengthInd | badPowerInd) = NaN;
    
end



end

