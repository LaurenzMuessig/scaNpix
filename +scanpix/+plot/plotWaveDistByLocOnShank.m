function plotWaveDistByLocOnShank(waveforms,coords)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


figSz  = [200 850];
screenSz = get(0,'ScreenSize');
hFig = figure('Units','pixels','Position',[0.1*screenSz(3) 0.1*screenSz(4) figSz]);
hAx = axes(hFig,'position',[0.25 0.05 0.7 0.9]);

% mean waveforms
meanWFs = cellfun(@(x) squeeze(nanmean(x,1)),waveforms,'uni',0);
[maxVals, maxInd] = cellfun(@(x) max(abs(min(x,[],1) - max(x,[],1))),meanWFs,'uni',0);
maxVal = max(cell2mat(maxVals));

% tempCh     = cell2mat(cellfun(@(x,y) x(y), app.mainWin.UserData.dataObj{currDsSelInd}.spikeData.spk_waveforms{currTrSelInd,2}(ind),maxInd,'uni',0));

nSamplesWaveform = size(meanWFs{1},1);
col   = scanpix.fxchange.cbrewer('qual', 'Set1', max([length(meanWFs), 3]), 'PCHIP' ); % cbrewer makes nice colormaps
xBase = [11 27 43 59; [0 1 2 3] .* (nSamplesWaveform/2)];

hold(hAx,'on')
for j = 1:length(tempWaves)
%     coords(j,:) = [app.mainWin.UserData.dataObj{currDsSelInd}.chanMap(currTrSelInd).xcoords(tempCh(j)) app.mainWin.UserData.dataObj{currDsSelInd}.cell_ID(tempCh(j),2)/15];
    
    xStart = xBase(2,coords(j,1)==xBase(1,:));
%     plot(hAx,(xStart:xStart+40)+((-5-5)*rand(1)+5),-((meanWFs{j}(:,maxInd{j})./maxVal)*5)+(coords(j,2)+((-3-3)*rand(1)+3)),'color',col(j,:),'LineWidth',2);
    plot(hAx,(xStart:xStart+nSamplesWaveform-1),-((meanWFs{j}(:,maxInd{j})./maxVal)*5)+coords(j,2),'color',col(j,:),'LineWidth',2);

end
hold(hAx,'off');
set(hAx,'xlim',[-3 107],'ylim',[min(coords(:,2))-5 max(coords(:,2))+5],'xtick',[],'ytick',ylim,'yticklabel',{'ventral';'dorsal'})  
 
end
