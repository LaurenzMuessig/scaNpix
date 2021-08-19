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

%%TO-DO:
% test on different resolution machines

% 
channelSpacing = 20;
scalingFact = (max(coords(:,2)) - min(coords(:,2)) ) / 15; % dividing by 15 seems to yield okay size for waveforms, irrespective of length of probe that is plotted 
figSz  = [400 850]; % hard coded

% scale axis if more than 600um of shank are plotted
shankDist2PLot = (max(coords(:,2)) - min(coords(:,2)));
if shankDist2PLot < 600
    yScaling = 1;
else
    yScaling = shankDist2PLot / 600; % note hard coded!
end
% open scrollable figure
screenSz = get(0,'ScreenSize');
hScroll              = scanpix.plot.createScrollPlot([0.1*screenSz(3) 0.1*screenSz(4) figSz]); 
hScroll.hFig.Name    = 'scaNpix::WaveformDistOnShank';
hAx                  = scanpix.plot.addAxisScrollPlot( hScroll,[0.15*hScroll.hFig.Position(3) 0.05*hScroll.hFig.Position(4) 0.7*hScroll.hFig.Position(3) round((figSz(2)-50) * yScaling)], [0.15*hScroll.hFig.Position(3) 0.05*hScroll.hFig.Position(4)]);

% need overall of max of waveform to scale data properly
[maxVals, maxInd] = cellfun(@(x) max(abs(min(x,[],1) - max(x,[],1))),meanWaveforms,'uni',0);
maxVals = max(cell2mat(maxVals));
            
nSamplesWaveform = size(meanWaveforms{1},1);
xBase  = [11 43 27 59; [0 2 1 3] .* (nSamplesWaveform/2)];

% if channel map is supplied we plot probe recording sites in background,
% otherwise just line as indication
if nargin < 3
    plot(hAx,[xBase(2,:)+nSamplesWaveform/2;xBase(2,:)+nSamplesWaveform/2],[min(coords(:,2))*ones(1,4);max(coords(:,2))*ones(1,4)],'k-');
else
    chanMap = sortrows(chanMap,2);
    startInd = find(chanMap(:,2) < min(coords(:,2)) ,1,'last' );
    endInd = find(chanMap(:,2) > max(coords(:,2)) ,1,'first' );
    tempChanMap = chanMap(startInd:endInd,:);
    [~,idx] = ismember(tempChanMap(:,1), xBase(1,:));
    plot(hAx,xBase(2,idx)+nSamplesWaveform/2,tempChanMap(:,2),'ks','MarkerFaceColor','k','MarkerSize',4);
end

col    = scanpix.fxchange.cbrewer('qual', 'Set1', max([length(meanWaveforms), 3]), 'PCHIP' ); % cbrewer makes nice colormaps
colInd = randperm(size(col,1)); % randomise so colours don't cluster in groups
% plot waveforms
hold(hAx,'on');
for j = 1:length(meanWaveforms)
   
    xStart = xBase(2,coords(j,1)==xBase(1,:)) + randi(round([-nSamplesWaveform/4 nSamplesWaveform/4]),1); % jitter by a quarter of n waveform samles in x 
    plot(hAx,(xStart:xStart+nSamplesWaveform-1),-((meanWaveforms{j}(:,maxInd{j})./maxVals)*scalingFact)+coords(j,2) + randi(round([-channelSpacing/2 channelSpacing/20]),1),'color',col(colInd(j),:),'LineWidth',2); % jitter by half of channel distance in y

end
hold(hAx,'off');

% set(hAx,'xlim',[-3 xBase(2,4)+nSamplesWaveform+3],'xtick',sort(xBase(2,:))+nSamplesWaveform/2,'xticklabel',sort(xBase(1,:)),'ylim',[min(coords(:,2))-scalingFact/4 max(coords(:,2))+scalingFact/2],'ytick',[]);%[min(coords(:,2)) max(coords(:,2))],'yticklabel',round([min(coords(:,2)),max(coords(:,2))]));  
set(hAx,'xlim',[-3 xBase(2,4)+nSamplesWaveform+3],'xtick',[],'ylim',[min(coords(:,2))-scalingFact/4 max(coords(:,2))+scalingFact/2],'ytick',[]);%[min(coords(:,2)) max(coords(:,2))],'yticklabel',round([min(coords(:,2)),max(coords(:,2))]));  

end
