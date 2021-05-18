function plotMapsGroup(maps,type,varargin)
%  
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% To DO
% add case where you might want to plot e.g. sac's along rate maps 
%%


%% parse input
defaultCellIDStr      = strrep(string(strcat('c_',num2str((1:size(maps{1}))'))),' ','');
defaultCMap           = 'jet';
defaultNSteps         = 11; 
defaultBaseWdthFigure = 185;      % pixel
defaultPlotSize       = [75 75];  % pixel
defaultPlotSep        = [40 30];  % pixel
defaultOffsetBase     = [60 50];  % pixel
defaultFigName        = 'MapGroup';  
% 
p = inputParser;
addOptional(p,'cellIDStr',defaultCellIDStr);
addParameter(p,'cmap',defaultCMap,@ischar);
addParameter(p,'nsteps',defaultNSteps,@isscalar);
addParameter(p,'basewidth',defaultBaseWdthFigure);
addParameter(p,'plotsize',defaultPlotSize);
addParameter(p,'plotsep',defaultPlotSep);
addParameter(p,'offsetbase',defaultOffsetBase);
addParameter(p,'figname',defaultFigName,@ischar);
parse(p,varargin{:});


%% plot
% set up figure 
screenSz       = get(0,'screensize');

% figure width
wdthScaleFact  = (length(maps)-1) * (p.Results.basewidth/1.5); % add this to fig width for any additional cell
figWidth       = p.Results.basewidth + wdthScaleFact;
% threshold size
if figWidth > 0.7*screenSz(3)
    figWidth = 0.7*screenSz(3);
end
figSize     = [0.1*screenSz(3) 0.1*screenSz(4) figWidth 0.7*screenSz(4)]; % 

% gather plot tiling
if length(maps) == 1 
    nPlotsPerRow = 6; % in case just 1 trial, we want to make a compact plot...
else
    nPlotsPerRow = length(maps); % ...if it's several trials, we plot nCells across trials
%     nPlotsPerRow = nPlotsPerRow + nPlotsPerRow*prms.plot_sAC;
end
% open plot
offsets   = p.Results.offsetbase;
hScroll = scanpix.plot.createScrollPlot( [figSize(1:2) nPlotsPerRow*p.Results.plotsize(1)+nPlotsPerRow*p.Results.plotsep(1)+3*offsets(1) figSize(4) ]  );
hScroll.hFig.Name = p.Results.figname;
% open a waitbar
nPlots    = length(maps)*length(maps{1});
plotCount = 1;
hWait     = waitbar(0); 

nRowPlots = 0;

for i = 1:length(maps{1})
    
    for j = 1:length(maps)
        
        waitbar(plotCount/nPlots,hWait,['Plotting ' type ' maps, just bare with me!']);
        
        % plot
        hAx = scanpix.plot.addAxisScrollPlot( hScroll, [offsets p.Results.plotsize], p.Results.plotsep );
        
        if strcmpi(type,'rate') || strcmpi(type,'spike') || strcmpi(type,'pos')
            scanpix.plot.plotRateMap(maps{j}{i},hAx,'colmap',p.Results.cmap,'nsteps',p.Results.nsteps)
        elseif  strcmp(type,'dir')
            scanpix.plot.plotDirMap(maps{j}{i},hAx);
        elseif  strcmpi(type,'sacs')
            imagesc(hAx,maps{j}{i}); colormap(hAx,jet);
            axis(hAx,'square');
            mapSz = size(maps{j}{i});
            set(hAx,'ydir','normal','xlim',[0 mapSz(2)],'ylim',[0 mapSz(1)]);
            axis(hAx,'off');
        end
        
        % plot peak rate
%         if ~strcmpi(type,'sacs')
            t = text(hAx);
            set(t,'Units','pixels','position',[8 -6],'String',sprintf('fRate=%.1f',nanmax(maps{j}{i}(:)) ),'FontSize',8 ); % harcoded text pos
%         end
        % plot cell ID string
        if j == 1
            t = text(hAx);
            set(t,'Units','pixels','position',[-40 hAx.Position(4)/2],'String',p.Results.cellIDStr{i},'FontSize',8,'Interpreter','none' ); % harcoded text pos
        end
        % update
        offsets(1) = offsets(1) + p.Results.plotsize(1) + 1.5*p.Results.plotsep(1);
        nRowPlots = nRowPlots + 1; %+ prms.plot_sAC;
        
        if nRowPlots >= nPlotsPerRow
            offsets(1) = p.Results.offsetbase(1);
            offsets(2) = offsets(2) + p.Results.plotsize(2) + p.Results.plotsep(2);
            nRowPlots = 0;
        end
%         offsets(2) = offsets(2) + p.Results.plotsize(2) + p.Results.plotsep(2);
        
        if i == length(maps{1})
            title(hAx,['Trial_' num2str(j)],'Interpreter','none');
        end
        plotCount = plotCount + 1;
    end
    % update and reset
%     offsets(1) = i * (p.Results.plotsize(1)+p.Results.plotsep(1)) + p.Results.offsetbase(1);
%     offsets(2) = p.Results.offsetbase(2);
end
close(hWait);

end

