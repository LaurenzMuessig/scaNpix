function plotWaveform(waveforms,varargin)
%
%% parse input
defaultMaxWaves         = 250;
defaulthAx              = [];
defaultPlotIndividWaves = true;
% 
p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'maxWaves',defaultMaxWaves,@isscalar);
addParameter(p,'plotIndWFs',defaultPlotIndividWaves,@islogical);
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
plot(hAx,meanWFs(:,maxInd),'color',[1 0 0],'linewidth',6);
hold(hAx, 'off');
set(hAx,'xlim',[1 size(meanWFs,1)],'ylim',axLims);

ylabel(hAx,'amplitude (\muV)');
    
end

