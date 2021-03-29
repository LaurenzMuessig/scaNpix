function plotWaveDistByLocOnShank(meanWaveforms,coords,chanMap)
% plotWaveDistByLocOnShank - plot the distribution of waveforms along the 
% shank of a neuropixel probe. All units have their amplitudes scaled to the
% max amplitude in data.
%
% Usage:
%       [templateDepths, maxChan] = scanpix.plot.plotWaveDistByLocOnShank(meanWaveforms,coords);
%       [templateDepths, maxChan] = scanpix.plot.plotWaveDistByLocOnShank(meanWaveforms,coords,chanMap);
%
% Inputs:   meanWaveforms   - cell array of mean waveforms to plot
%           coords          - location for units on shank as x,y pairs (for y it's best to use COM rather than max channel to avoid quantisation of data)
%           chanMap         - (optional), full channel map as x,y pairs
%
% Outputs:  
%
% LM 2021

% hard coded at the moment - should make it scale on different screens
scalingFact = 100;
figSz  = [200 850];
screenSz = get(0,'ScreenSize');
hFig = figure('Units','pixels','Position',[0.1*screenSz(3) 0.1*screenSz(4) figSz]);
hAx = axes(hFig,'position',[0.25 0.05 0.7 0.9]);

% need overall of max of waveform to scale data properly
[maxVals, maxInd] = cellfun(@(x) max(abs(min(x,[],1) - max(x,[],1))),meanWaveforms,'uni',0);
maxVals = max(cell2mat(maxVals));
            
nSamplesWaveform = size(meanWaveforms{1},1);
xBase  = [11 43 27 59; [0 2 1 3] .* (nSamplesWaveform/2)];

% if channel map is supplied we plot probe recording sites in background,
% otherwise just line as indication
if nargin == 3
    plot(hAx,[xBase(2,:)+nSamplesWaveform/2;xBase(2,:)+nSamplesWaveform/2],[min(coords(:,2))*ones(1,4);max(coords(:,2))*ones(1,4)],'k:');
else
    chanMap = sortrows(chanMap,2);
    startInd = find(chanMap(:,2) < min(coords(:,2)) ,1,'last' );
    endInd = find(chanMap(:,2) > max(coords(:,2)) ,1,'first' );
    tempChanMap = chanMap(startInd:endInd,:);
    [~,idx] = ismember(tempChanMap(:,1), xBase(1,:));
    plot(hAx,xBase(2,idx)+nSamplesWaveform/2,tempChanMap(:,2),'ks','MarkerFaceColor','k');
end

col    = scanpix.fxchange.cbrewer('qual', 'Set1', max([length(meanWaveforms), 3]), 'PCHIP' ); % cbrewer makes nice colormaps
colInd = randperm(size(col,1)); % randomise so colours don't cluster in groups
% plot waveforms
hold(hAx,'on');
for j = 1:length(meanWaveforms)
   
    xStart = xBase(2,coords(j,1)==xBase(1,:));
    plot(hAx,(xStart:xStart+nSamplesWaveform-1),-((meanWaveforms{j}(:,maxInd{j})./maxVals)*scalingFact)+coords(j,2),'color',col(colInd(j),:),'LineWidth',2);

end
hold(hAx,'off');

set(hAx,'xlim',[-3 xBase(2,4)+nSamplesWaveform+3],'xtick',sort(xBase(2,:))+nSamplesWaveform/2,'xticklabel',sort(xBase(1,:)),'ylim',[min(coords(:,2))-scalingFact/4 max(coords(:,2))+scalingFact],'ytick',[min(coords(:,2)) max(coords(:,2))],'yticklabel',round([min(coords(:,2)),max(coords(:,2))]));  
end
