function [peakFreq, maxPower, s2n, bestEEGInd] = lfpPowerSpec(eeg, freqBand, sampleRate, options)
% eegPowerSpec - Calculate power spectrum for eeg using the FFT 
% Extract peak frequency, power and signal-to-noise ratios in specific band
%
% Usage:
%
%       [peakFreq, maxPower, s2n, bestEEGInd] = eegPowerSpec(eeg, freqBand, varargin)
%       [peakFreq, maxPower, s2n, bestEEGInd] = eegPowerSpec(eeg, freqBand, optionalInputStruct );
%       [peakFreq, maxPower, s2n, bestEEGInd] = eegPowerSpec(eeg, freqBand, 'inputName', inputVal, .. etc .. );
%
% Inputs:   eeg      - nSamples x nEEG array of eeg voltages 
%           freqBand - [lowerBound upperBound] of frequency band from which we want to get frequency with max power (in Hz)
%
%
% Options:
%
%   'eegFs',          250,  - sample rate for EEGs
%   'maxFreq',        [],   - max frequency for spectrogram in Hz - if empty we'll use Nyquist.
%   'smWinSz',        2     - size of smoothing window for spectrogram (in Hz)
%   's2nBand',        0.5   - freq band around peak freq for signal2noise ratio calculation
%   'plotSpec',       0     - flag if you want to plot the spectrogram
%
% Outputs:  peakFreq   - nEEG x 1 frequencies @ max. power in [freqBand] 
%           maxPower   - nEEG x 1 max. power in [freqBand] 
%           s2n        - nEEG x 1 signal-to-noise ratios in [freqBand] 
%           bestEEGInd - index of best EEG ( i.e. max(s2n) ) 
%
% This is inspired by 'eeg_powerspec' in Scan (by TW), but will work on an array
% of EEGs from a given trial instead of just one at a time
%
% LM 2020

%%
arguments
    eeg {mustBeNumeric}
    freqBand (1,2) {mustBeNumeric}
    sampleRate  (1,1) {mustBeNumeric}
    options.maxFreq {mustBeNumeric} = [];
    options.smWinSz (1,1) {mustBeNumeric} = 2;
    options.s2nBand (1,1) {mustBeNumeric} = 0.5;
    options.plotSpec (1,1) {mustBeNumericOrLogical} = false;
end

%%
% a bit more input checking...
% make sure array is nSamples x nEEGs 
if size(eeg,1) < size(eeg,2)
    eeg = eeg';
end
%
if isempty(options.maxFreq)
    options.maxFreq = sampleRate/2; % Nyquist freq
end

%% remove NaNs
eeg       = eeg( ~isnan(eeg(:,1)), : ); % we are assuming that NaNs will be the same in every eeg after e.g. speed filtering 

%% remove DC shift (i.e. detrend)
eeg       = bsxfun(@minus,eeg,mean(eeg,1,'omitnan'));

%% Fourier transform 
nFFT      = nextpow2(length(eeg));
fourTrans = fft(eeg, 2^nFFT, 1);
% Frequencies at each power spectrum sample point
freq     = ( 0 : length(fourTrans) - 1) * sampleRate / length(fourTrans);
freq     = freq( freq <= options.maxFreq ); % drop redundant frequencies

%% power
tempPower = abs(fourTrans).^2 / length(fourTrans);  % same as power = fourTrans .* conj(fourTrans) / length(fourTrans), but for matrix abs(fourTrans).^2 is faster 
tempPower = tempPower(1:length(freq), : ); % drop power for redundant frequencies
% smooth power
kernSz    = find(freq <= options.smWinSz, 1 ,'last');
kernel    = fspecial( 'Gaussian', [kernSz 1], kernSz/3 ) ; % TW uses sigma = kernSz/3
% smooth - same as imfilter(tempPower(:,i), kernel, 'replicate'), but imfilter gets extremely slow for long
% recordings (kernel size in bins scales with nFFT; e.g. ~18min vs 0.16s with conv2 for a 15min trial @ 2.5kHz)
kernCtr   = floor((kernSz+1)/2); % same kernel centre as imfilter
padded    = [repmat(tempPower(1,:),kernCtr-1,1); tempPower; repmat(tempPower(end,:),kernSz-kernCtr,1)]; % replicate edges
power     = conv2(padded, kernel, 'valid'); % smooths all EEGs (columns) at once

%% find peak freq
freqBandPower       = power( freq > freqBand(1) & freq < freqBand(2), : );
freqBandFreqs       = freq( freq > freqBand(1) & freq < freqBand(2) );
[maxPower, freqInd] = max( freqBandPower, [], 1 );
peakFreq            = freqBandFreqs(freqInd);
% if there isn't a clear (e.g. theta) peak in spectrum the peak freq will typically be the first index of band (lower F = higher power)
% in that case we want to find the local max in frequency band with highest power
if any(freqInd == 1)
    temp_freqBandPower                   = freqBandPower( :, freqInd == 1 );
    localMaxInd                          = islocalmax(temp_freqBandPower, 1); % local maxima
    localMaxInd(1, ~any(localMaxInd, 1)) = true; % protect against cases w/o local max - restore original value
    temp_freqBandPower(~localMaxInd)     = 0;
    [tempMaxPow, tempInd]                = max(temp_freqBandPower, [], 1);
    tempPeakFreq                         = freqBandFreqs(tempInd);
    % adjust max power and peak freq
    maxPower(freqInd == 1)               = tempMaxPow;
    peakFreq(freqInd == 1)               = tempPeakFreq;
end


%% get s2n ratio
peakFreqInd = freq > peakFreq' - options.s2nBand & freq < peakFreq' + options.s2nBand;
s2n         = mean( reshape( power(peakFreqInd'), [], size(power,2) ), 1, 'omitnan' ) ./ mean( reshape( power(~peakFreqInd'), [], size(power,2) ), 1, 'omitnan' ); % might make sense to exclude all f<1-2Hz for this?

%% best EEG index
[~, bestEEGInd] = max(s2n);

%% plot
if options.plotSpec
    figure;
    otherInd = setdiff(1:size(power,2), bestEEGInd);
    if ~isempty(otherInd) % no other EEGs if only 1 EEG supplied
        plot(freq, power(:,otherInd)', 'k-');
    end
    hold on
    plot(freq,power(:,bestEEGInd)','r-','linewidth', 3);
    line([peakFreq(bestEEGInd)+0.1*options.maxFreq peakFreq(bestEEGInd)],[maxPower(bestEEGInd)*1.1 maxPower(bestEEGInd)],'color',[0.5 0.5 0.5]);
    hold off
    set(gca, 'ylim', [0 maxPower(bestEEGInd)*1.3], 'xlim', [0 options.maxFreq]);
    text(peakFreq(bestEEGInd)+0.1*options.maxFreq,maxPower(bestEEGInd)*1.1,sprintf('Freq=%.2f',peakFreq(bestEEGInd)));
    xlabel('Frequency (Hz)');
    ylabel('Power (\muV^2)');
end

end


