function plotWaveform(waveforms,options)
% plotWaveform - plot wave form of a cell. ATM will only plot the waveform
% on max channel
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotWaveform( waveforms ) 
%           scanpix.plot.plotWaveform( waveforms, hAx )
%           scanpix.plot.plotWaveform( __ ,'name',value,.... )
%
%  Inputs:  
%           waveforms - array of raw waveform samples (nSamples-by-nSamplesPerWaveform-by-nChannels)
%           varargin  - optional inputs
%                     - axis handle and/or
%                     - Name-Value pairs for 'maxWaves' (plot this many waveform samples)
%                       and/or 'plotIndWFs' (flag to plot individual waveforms)
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% parse input
arguments
    waveforms {mustBeNumeric} 
    options.ax {ishghandle(options.ax, 'axes')}
    options.maxWaves (1,1) {mustBeNumeric} = 250; 
    options.linew (1,1) {mustBeNumeric} = 3; 
    options.plotAllWaves (1,1) {mustBeNumericOrLogical} = true;
end

% 
if isempty(waveforms)
    return;
end
%
if ~isfield(options,'ax')
    options.ax = axes;
end


%% plot
if size(waveforms,3) == 1
    % in case of 1 spike, duplicate that to make code compatible 
%     tempWF = cat(3,waveforms,waveforms); % TEst this does the right thing for 1 spike case
%     if size(tempWF,3) ~= 2
%         waveforms = shiftdim(tempWF,2);
%     end
%     waveforms = shiftdim(cat(3,waveforms,waveforms),2);
end


meanWFs = squeeze(mean(waveforms,1,'omitnan'));
if size(meanWFs,1) == 1; meanWFs = meanWFs'; end

[~, maxInd] = max(abs(min(meanWFs,[],1) - max(meanWFs,[],1)));
if options.plotAllWaves
    axLims  = [min(min(waveforms(:,:,maxInd))) max(max(waveforms(:,:,maxInd)))];
else
    axLims  = [min(meanWFs(:)) max(meanWFs(:))];
end

nSamples = size(waveforms,2);
if size(waveforms,1) > options.maxWaves
    selInd    = round(linspace(1,size(waveforms,1),options.maxWaves));
    waveforms = [squeeze(waveforms(selInd,:,maxInd)) nan(options.maxWaves,1)]; % adding NaNs will results in plotting a sngle line object
    nWaves    = options.maxWaves;
else
    nWaves    = size(waveforms,1);
    waveforms = [squeeze(waveforms(:,:,maxInd)) nan(nWaves,1)];
end
waveforms     = waveforms';

if options.plotAllWaves; plot(options.ax,repmat([1:nSamples,NaN],1,nWaves)',waveforms(:),'color',[0.5 0.5 0.5],'linewidth',.3); end
hold(options.ax, 'on');
plot(options.ax,meanWFs(:,maxInd),'color',[1 0 0],'linewidth',options.linew);
hold(options.ax, 'off');
set(options.ax,'xlim',[1 size(meanWFs,1)],'ylim',axLims);

ylabel(options.ax,'amplitude (\muV)');
    
end

