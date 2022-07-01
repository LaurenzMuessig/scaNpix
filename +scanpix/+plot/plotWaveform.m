function plotWaveform(waveforms,varargin)
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
defaultMaxWaves         = 250;
defaulthAx              = [];
defaultPlotIndividWaves = true;
defaultLineWidth        = 3;
% 
p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'maxWaves',defaultMaxWaves,@isscalar);
addParameter(p,'plotIndWFs',defaultPlotIndividWaves,@islogical);
addParameter(p,'linew',defaultLineWidth);
parse(p,varargin{:});
% 
if isempty(p.Results.ax)
    hAx = axes;
else
    hAx = p.Results.ax;
end

if isempty(waveforms)
    return;
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


meanWFs = squeeze(nanmean(waveforms,1));
if size(meanWFs,1) == 1; meanWFs = meanWFs'; end

[~, maxInd] = max(abs(min(meanWFs,[],1) - max(meanWFs,[],1)));
if p.Results.plotIndWFs
    axLims  = [min(min(waveforms(:,:,maxInd))) max(max(waveforms(:,:,maxInd)))];
else
    axLims  = [min(meanWFs(:)) max(meanWFs(:))];
end

nSamples = size(waveforms,2);
if size(waveforms,1) > p.Results.maxWaves
    selInd = round(linspace(1,size(waveforms,1),p.Results.maxWaves));
    waveforms = [squeeze(waveforms(selInd,:,maxInd)) nan(p.Results.maxWaves,1)]; % adding NaNs will results in plotting a sngle line object
    nWaves = p.Results.maxWaves;
else
    nWaves = size(waveforms,1);
    waveforms = [squeeze(waveforms(:,:,maxInd)) nan(nWaves,1)];
end
waveforms = waveforms';

if p.Results.plotIndWFs; plot(hAx,repmat([1:nSamples,NaN],1,nWaves)',waveforms(:),'color',[0.5 0.5 0.5],'linewidth',.3); end
hold(hAx, 'on');
plot(hAx,meanWFs(:,maxInd),'color',[1 0 0],'linewidth',p.Results.linew);
hold(hAx, 'off');
set(hAx,'xlim',[1 size(meanWFs,1)],'ylim',axLims);

ylabel(hAx,'amplitude (\muV)');
    
end

