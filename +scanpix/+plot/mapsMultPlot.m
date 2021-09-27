function mapsMultPlot(data,type,varargin)
%  
% This is a generic plotting routine to make a scrollable plot for any
% arbitrary combination of maps. We make the following assumption that 
% data = cell(1,nTrials,nMapTypes) 
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 

%% parse input
defaultCellIDStr      = strrep(string(strcat('c_',num2str((1:size(data{1}))'))),' ','');
defaultCMap           = 'jet';
defaultNSteps         = 11; 
defaultBaseWdthFigure = 185;      % pixel
defaultPlotSize       = [75 75];  % pixel
defaultPlotSep        = [40 30];  % pixel
defaultOffsetBase     = [60 50];  % pixel
defaultFigName        = 'multPlot';  
% 
p = inputParser;
addOptional(p,'cellIDStr',defaultCellIDStr,@isstring);
addParameter(p,'cmap',defaultCMap,@ischar);
addParameter(p,'nsteps',defaultNSteps,@isscalar);
addParameter(p,'basewidth',defaultBaseWdthFigure);
addParameter(p,'plotsize',defaultPlotSize);
addParameter(p,'plotsep',defaultPlotSep);
addParameter(p,'offsetbase',defaultOffsetBase);
addParameter(p,'figname',defaultFigName,@ischar);
parse(p,varargin{:});

% some sanity checcks should go here

%% plot
% set up figure 
screenSz       = get(0,'screensize');

% figure width
wdthScaleFact  = (numel(data)-1) * (p.Results.basewidth/1.5); % add this to fig width for any additional trial and plot type
figWidth       = p.Results.basewidth + wdthScaleFact;
% threshold size
if figWidth > 0.7*screenSz(3)
    figWidth = 0.7*screenSz(3);
end
figSize     = [0.1*screenSz(3) 0.1*screenSz(4) figWidth 0.7*screenSz(4)]; % 

% gather plot tiling
if numel(data) == 1 
    nPlotsPerRow = 6; % in case just 1 trial, we want to make a compact plot...
    noGridMode   = false;
else
    nPlotsPerRow = numel(data); % ...if it's several trials/plot types, we plot nCells across those
    noGridMode   = true;
end
% open plot
offsets              = p.Results.offsetbase;
hScroll              = scanpix.plot.createScrollPlot( [figSize(1:2) nPlotsPerRow*p.Results.plotsize(1)+nPlotsPerRow*p.Results.plotsep(1)+3*offsets(1) figSize(4) ]  );
hScroll.hFig.Name    = p.Results.figname;
hScroll.hFig.Visible = 'off'; % hiding figure speeds up plotting by a fair amount
% open a waitbar
nPlots    = prod([size(data) length(data{1})]);
plotCount = 1;
hWait     = waitbar(0); 

nRowPlots = 0;

for i = 1:length(data{1})
    
    for j = 1:size(data,3)
        
        plotPeakRateFlag = false;
        
        for k = 1:size(data,2)
            
            waitbar(plotCount/nPlots,hWait,'Making your precious figure, just bare with me!');
            
            % plot
            hAx = scanpix.plot.addAxisScrollPlot( hScroll, [offsets p.Results.plotsize], p.Results.plotsep );
            
            if strcmpi(type{j},'rate') || strcmpi(type{j},'spike') || strcmpi(type{j},'pos')
                scanpix.plot.plotRateMap(data{1,k,j}{i},hAx,'colmap',p.Results.cmap,'nsteps',p.Results.nsteps)
                plotPeakRateFlag = true;
            elseif  strcmp(type{j},'dir')
                scanpix.plot.plotDirMap(data{1,k,j}{i},hAx);
                plotPeakRateFlag = true;
            elseif  strcmp(type{j},'lin')
                 scanpix.plot.plotLinMaps( data{1,k,j}{i}, hAx);
            elseif  strcmpi(type{j},'sacs')
                mapSz = size(data{1,k,j}{i});
                imagesc(hAx,'CData',data{k}{i}); colormap(hAx,jet);
                axis(hAx,'square');
                
                set(hAx,'xlim',[0 mapSz(2)],'ylim',[0 mapSz(1)]);
                axis(hAx,'off');
            elseif  strcmpi(type{j},'speed')
                scanpix.plot.plotSpeedMap(data{1,k,j}{i},hAx);
            elseif  strcmpi(type{j},'custom')
                %%% PROBABLY WOULD REQUIRE SUPPLYING ANON FNCT TO PLOT
            end
            
            % plot peak rate
            if plotPeakRateFlag
                t = text(hAx);
                set(t,'Units','pixels','position',[8 -6],'String',sprintf('peakFR=%.1f',nanmax(data{1,k,j}{i}(:)) ),'FontSize',8 ); % harcoded text pos
            end
            % plot cell ID string
            if j == 1 && k == 1
                t = text(hAx);
                set(t,'Units','pixels','position',[-40 hAx.Position(4)/2],'String',p.Results.cellIDStr{i},'FontSize',8,'Interpreter','none' ); % harcoded text pos
            end
            % update
            offsets(1) = offsets(1) + p.Results.plotsize(1) + 1.5*p.Results.plotsep(1);
            nRowPlots = nRowPlots + 1; %;
            
            if nRowPlots >= nPlotsPerRow
                offsets(1) = p.Results.offsetbase(1);
                offsets(2) = offsets(2) + p.Results.plotsize(2) + p.Results.plotsep(2);
                nRowPlots = 0;
            end
            
            if i == length(data{1}) && noGridMode
                title(hAx,['Trial_' num2str(k) '-' type{j}],'Interpreter','none');
            end
            plotCount = plotCount + 1;
            %disp(plotCount);
        end
    end
end

close(hWait);

hScroll.hFig.Visible = 'on';
end

