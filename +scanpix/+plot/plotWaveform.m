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
    waveforms = shiftdim(cat(3,waveforms,waveforms),2);
end

meanWFs = squeeze(nanmean(waveforms,1));
if size(meanWFs,1) == 1; meanWFs = meanWFs'; end

[~, maxInd] = max(abs(min(meanWFs,[],1) - max(meanWFs,[],1)));
axLims  = [min(min(waveforms(:,:,maxInd))) max(max(waveforms(:,:,maxInd)))];
if size(waveforms,1) > p.Results.maxWaves
    selInd = round(linspace(1,size(waveforms,1),p.Results.maxWaves));
    waveforms = squeeze(waveforms(selInd,:,maxInd));
else
    waveforms = squeeze(waveforms(:,:,maxInd));
end


if p.Results.plotIndWFs; plot(hAx,waveforms','color',[0.5 0.5 0.5],'linewidth',.3); end
hold(hAx, 'on');
plot(hAx,meanWFs(:,maxInd),'color',[1 0 0],'linewidth',p.Results.linew);
hold(hAx, 'off');
set(hAx,'xlim',[1 size(meanWFs,1)],'ylim',axLims);

ylabel(hAx,'amplitude (\muV)');
    
end

