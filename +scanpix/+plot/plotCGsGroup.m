function plotCGsGroup(spikeTimes,binSz,lag,trialDur,varargin)
% plot all possible temporal cross-correlograms (CGs) based on a cell array  
% with spike times of different cells. Useful to check if clusters might
% have to be merged. We'll make a figure that is scrollable.
% Note: Can be slow if you try to plot too many cells. 
%
% Usage:    
%           scanpix.plot.plotCGsGroup(spikeTimes,binSz,lag,trialDur)
%           scanpix.plot.plotCGsGroup(spikeTimes,binSz,lag,trialDur,cell_IDs)
%
% Inputs:   
%           spikeTimes - cell array of spike times of different single units
%           binSz      - bin size for correlograms in seconds
%           lag        - max lag for correlograms in seconds
%           trialDur   - trial duration in seconds
%           cell_IDs   - cell array of cell ID strings (optional)      
%
% Outputs:  
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TO DO
% text positioning doesn't work well across different screens


%% params
rMaps          = {};
cell_IDs      = strcat({'cell_'},num2str((1:length(spikeTimes))'));
figName       = 'scaNpix::CG_overview';
plotSize      = [75 75];
plotSep       = [15 20];
offsetBase    = [50 40];
% baseSzFigure  = [180 170];
plotRMaps     = false;
groupInd      = [];

p = inputParser;
addOptional(p,'maps',rMaps,@iscell);
addParameter(p,'cellIDStr',cell_IDs,@(x) isstring(x) || iscell(x));
addParameter(p,'figname',figName,@ischar);
addParameter(p,'plotsize',plotSize);
addParameter(p,'plotsep',plotSep);
addParameter(p,'offsetbase',offsetBase);
% addParameter(p,'baseSz',baseSzFigure);
addParameter(p,'plotmaps',plotRMaps,@islogical);
addParameter(p,'groupInd',groupInd,@(x) islogical(x) || isempty(x));
parse(p,varargin{:});


if p.Results.plotmaps
    if isempty(p.Results.maps); error('You need to supply maps if you want to plot them. Not sure why I have to tell you that...'); end
%     offsetBase   = p.Results.offsetbase +  p.Results.plotsize + p.Results.plotsep;
else
%     offsetBase   = p.Results.offsetbase ;
end


%% plot
% screenSz      = get(0,'screensize');
% wdthScaleFact = (length(spikeTimes)-1) * (p.Results.baseSz(1)/2); % add this to fig width for any additional cell
% figWdth       = p.Results.baseSz(1) + wdthScaleFact;
% % width
% if figWdth > 0.7*screenSz(3)
%     figWdth = 0.7*screenSz(3);
% end
% hghtScaleFact = (length(spikeTimes)-1) * (p.Results.baseSz(2)/2); % add this to fig width for any additional cell
% figHght       = p.Results.baseSz(2) + hghtScaleFact;
% % height
% if figHght > 0.7*screenSz(4)
%     figHght = 0.7*screenSz(4);
% end
% % final size
% figSize     = [0.1*screenSz(3) 0.1*screenSz(4) figWdth figHght ]; 
% % open scrollable figure
% hScroll              = scanpix.plot.createScrollPlot(figSize); 
% hScroll.hFig.Name    = p.Results.figname;
% hScroll.hFig.Visible = 'off'; % hiding figure speeds up plotting by a fair amount
% %wait bar
% hWait         = waitbar(0); 
% nPlots        = (length(spikeTimes)*(length(spikeTimes)+1))/2 + p.Results.plotmaps * length(spikeTimes) * 2;
% plotCount     = 1;
% %
% offsets   = offsetBase;

%wait bar
hWait         = waitbar(0); 
nPlots        = (length(spikeTimes)*(length(spikeTimes)+1))/2 + p.Results.plotmaps * length(spikeTimes) * 2;
plotCount     = 1;

[axArr, hScroll] = scanpix.plot.multPlot([length(spikeTimes) length(spikeTimes)],'plotsize',p.Results.plotsize,'plotsep',p.Results.plotsep,'offsetbase',p.Results.offsetbase,'figname',p.Results.figname);
% hScroll.hFig.Visible = 'off';

for i = 1:length(spikeTimes)
    % wait bar
    waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
    % 
    if p.Results.plotmaps
        hAx        = scanpix.plot.addAxisScrollPlot( hScroll,[offsets+[0 p.Results.plotsep(2)]+[0 p.Results.plotsize(2)] p.Results.plotsize], p.Results.plotsep);
        scanpix.plot.plotRateMap(p.Results.maps{i},hAx);
        text(hAx,-12,0.45*max(get(hAx,'ylim')),p.Results.cellIDStr{i},'Interpreter','none');
    end
    
    % plot AC
%     hAx       = scanpix.plot.addAxisScrollPlot( hScroll,[offsets p.Results.plotsize], p.Results.plotsep);
    scanpix.analysis.spk_crosscorr(spikeTimes{i},'AC',binSz,lag,trialDur,'plot',axArr{i,1}); % autocorr
    yAxlim    = get(axArr{i,1},'ylim');
    if i ~= 1
        set(axArr{i,1},'xtick',[-lag 0 lag],'xticklabel',{''});
    else
        set(axArr{i,1},'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
    end
    % update offsets
%     offsets(1) = offsets(1) + p.Results.plotsize(1) + p.Results.plotsep(1);
    % plot headers
    if i==1 || ~p.Results.plotmaps 
        text(axArr{i,1},-lag-0.15,0.5*max(yAxlim),p.Results.cellIDStr{i},'Interpreter','none');
    end
    text(axArr{i,1},-0.3*lag,1.15*max(yAxlim),p.Results.cellIDStr{i},'Interpreter','none');
    plotCount = plotCount + 1;
    for j = i+1:length(spikeTimes)
        waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
        % add axis
%         hAx    = scanpix.plot.addAxisScrollPlot( hScroll,[offsets p.Results.plotsize], p.Results.plotsep );    
        % plot cross corr
        scanpix.analysis.spk_crosscorr(spikeTimes{i},spikeTimes{j},binSz,lag,trialDur,'plot',axArr{i,j}); % crosscorr
        yAxlim = get(axArr{i,j},'ylim');
         % plot comparison ID
        text(axArr{i,j},-lag-0.0025,1.15*max(yAxlim),[p.Results.cellIDStr{i} ' v ' p.Results.cellIDStr{j}],'Interpreter','none','color','r');
        if i~=1
            set(axArr{i,j},'xtick',[-lag 0 lag],'xticklabel',{''});
        else
            set(axArr{i,j},'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
        end
        if j == length(spikeTimes)
            if p.Results.plotmaps
                offsets(1) = offsets(1) + p.Results.plotsize(1) + p.Results.plotsep(1);
                hAx        = scanpix.plot.addAxisScrollPlot( hScroll,[offsets p.Results.plotsize], p.Results.plotsep);
                scanpix.plot.plotRateMap(p.Results.maps{i},hAx);
                text(hAx,max(get(hAx,'xlim'))+3,0.45*max(get(hAx,'ylim')),p.Results.cellIDStr{i},'Interpreter','none');
                plotCount  = plotCount + 1;
            else
                text(axArr{i,j},1.1*lag,0.5*max(yAxlim),p.Results.cellIDStr{i},'Interpreter','none');
            end
        end
        % update offsets    
%         offsets(1) = offsets(1) + p.Results.plotsize(1) + p.Results.plotsep(1);
        plotCount  = plotCount + 1;
    end
%     % update offsets  
%     offsets(2)     = offsets(2) + p.Results.plotsize(2) + p.Results.plotsep(2);
%     offsets(1)     = i * (p.Results.plotsize(1) + p.Results.plotsep(1)) + offsetBase(1);
end

if p.Results.plotmaps
    hAx        = scanpix.plot.addAxisScrollPlot( hScroll,[offsets - [0 p.Results.plotsize(2)] - [0 p.Results.plotsep(2)], p.Results.plotsize] , p.Results.plotsep);
    scanpix.plot.plotRateMap(p.Results.maps{i},hAx);
    text(hAx,max(get(hAx,'xlim'))+3,0.45*max(get(hAx,'ylim')),p.Results.cellIDStr{i},'Interpreter','none');
else
    text(axArr{i,j},1.1*lag,0.5*max(yAxlim),p.Results.cellIDStr{i},'Interpreter','none');
end


close(hWait);
% hScroll.hFig.Visible = 'on';
end


