function plotMapsGroup(maps,type,varargin)
%  
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% To DO
% ADD PLOT AS GRID IF ONLY 1 TRIAL
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
% open plot
hScroll     = scanpix.plot.createScrollPlot(figSize); % open scrollable figure
hScroll.hFig.Name = p.Results.figname;

offsets   = p.Results.offsetbase;
% open a waitbar
nPlots    = length(maps)*length(maps{1});
plotCount = 1;
hWait     = waitbar(0); 

for i = 1:length(maps)
    
    for j = 1:length(maps{i})
        
        waitbar(plotCount/nPlots,hWait,'Plotting rate maps, just bare with me!');
        
        % plot
        hAx = scanpix.plot.addAxisScrollPlot( hScroll, [offsets p.Results.plotsize], p.Results.plotsep );
        
        if strcmp(type,'rate') || strcmp(type,'spike') || strcmp(type,'pos')
            scanpix.maps.plotRateMap(maps{i}{j},hAx,'colmap',p.Results.cmap,'nsteps',p.Results.nsteps)
        elseif  strcmp(mapType,'dir')
            scanpix.maps.plotDirMap(maps{i}{j},hAx);
        end
        
        % plot peak rate
        t = text(hAx);
        set(t,'Units','pixels','position',[8 -6],'String',sprintf('fRate=%.1f',nanmax(maps{i}{j}(:)) ),'FontSize',8 ); % harcoded text pos
        % plot cell ID string
        if i == 1
            t = text(hAx);
            set(t,'Units','pixels','position',[-40 hAx.Position(4)/2],'String',p.Results.cellIDStr{j},'FontSize',8,'Interpreter','none' ); % harcoded text pos
        end
        % update
        offsets(2) = offsets(2) + p.Results.plotsize(2) + p.Results.plotsep(2);
        
        if j == length(maps{i})
            title(hAx,['Trial_' num2str(i)],'Interpreter','none');
        end
        plotCount = plotCount + 1;
    end
    % update and reset
    offsets(1) = i * (p.Results.plotsize(1)+p.Results.plotsep(1)) + p.Results.offsetbase(1);
    offsets(2) = p.Results.offsetbase(2);
end
close(hWait);

end

