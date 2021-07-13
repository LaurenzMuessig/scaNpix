function plotCGsGroup(spikeTimes,binSz,lag,trialDur,cell_IDs)
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


% hard coded properties
plotSize       = [75 75];
plotSep        = [15 20];
offsetBase     = [50 40];
baseWdthFigure = 180;
baseHghtFigure = 170;

%% params
if nargin < 5
    cell_IDs = strcat({'cell_'},num2str((1:length(spikeTimes))'));
end

%% plot
screenSz      = get(0,'screensize');
wdthScaleFact = (length(spikeTimes)-1) * (baseWdthFigure/2); % add this to fig width for any additional cell
figWdth       = baseWdthFigure + wdthScaleFact;
% width
if figWdth > 0.7*screenSz(3)
    figWdth = 0.7*screenSz(3);
end
hghtScaleFact = (length(spikeTimes)-1) * (baseHghtFigure/2); % add this to fig width for any additional cell
figHght       = baseHghtFigure + hghtScaleFact;
% height
if figHght > 0.7*screenSz(4)
    figHght = 0.7*screenSz(4);
end
% final size
figSize     = [0.1*screenSz(3) 0.1*screenSz(4) figWdth figHght ]; 
% open scrollable figure
hScroll       = scanpix.plot.createScrollPlot(figSize); 
hScroll.hFig.Visible = 'off'; % hiding figure speeds up plotting by a fair amount
%wait bar
hWait         = waitbar(0); 
nPlots        = (length(spikeTimes)*(length(spikeTimes)+1))/2;
plotCount     = 1;
offsets       = offsetBase;
for i = 1:length(spikeTimes)
    % wait bar
    waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
    % plot AC
    hAx       = scanpix.plot.addAxisScrollPlot( hScroll,[offsets plotSize], plotSep);
    scanpix.analysis.spk_crosscorr(spikeTimes{i},'AC',binSz,lag,trialDur,'plot',hAx); % autocorr
    yAxlim    = get(hAx,'ylim');
    if i ~= 1
        set(hAx,'xtick',[-lag 0 lag],'xticklabel',{''});
    else
        set(hAx,'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
    end
    % update offsets
    offsets(1) = offsets(1) + plotSize(1) + plotSep(1);
    % plot headers
    text(hAx,-lag-0.025,0.5*max(yAxlim),cell_IDs{i},'Interpreter','none');
    text(hAx,-0.3*lag,1.15*max(yAxlim),cell_IDs{i},'Interpreter','none');
    plotCount = plotCount + 1;
    for j = i+1:length(spikeTimes)
        waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
        % add axis
        hAx    = scanpix.plot.addAxisScrollPlot( hScroll,[offsets plotSize], plotSep );    
        % plot cross corr
        scanpix.analysis.spk_crosscorr(spikeTimes{i},spikeTimes{j},binSz,lag,trialDur,'plot',hAx); % crosscorr
        yAxlim = get(hAx,'ylim');
         % plot comparison ID
        text(hAx,-lag-0.0025,1.15*max(yAxlim),[cell_IDs{i} ' v ' cell_IDs{j}],'Interpreter','none','color','r');
        if i~=1
            set(hAx,'xtick',[-lag 0 lag],'xticklabel',{''});
        else
            set(hAx,'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
        end
        if j == length(spikeTimes)
            text(hAx,1.1*lag,0.5*max(yAxlim),cell_IDs{i},'Interpreter','none');   
        end
        % update offsets    
        offsets(1) = offsets(1) + plotSize(1) + plotSep(1);
        plotCount  = plotCount + 1;
    end
    % update offsets  
    offsets(2)     = offsets(2) + plotSize(2) + plotSep(2);
    offsets(1)     = i * (plotSize(1) + plotSep(1)) + offsetBase(1);
end
close(hWait);
hScroll.hFig.Visible = 'on';
end


