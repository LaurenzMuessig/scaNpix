function mapsMultPlot(data,type,varargin)
% This is a generic plotting routine to make a scrollable plot for any
% arbitrary combination of maps. We make the following assumption that 
% data = cell(1,nTrials,nMapTypes) 
% This fnct is a specialised version of scanpix.plot.multPlot  
% package: scanpix.plot
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 

%% parse input
defaultCellIDStr      = strrep(string(strcat('c_',num2str((1:size(data{1}))'))),' ','');
defaultCMap           = 'jet';
defaultNSteps         = 11; 
defaultPlotSize       = [75 75];  % pixel
defaultPlotSep        = [30 30];  % pixel
defaultOffsetBase     = [60 50];  % pixel
defaultFigName        = 'multPlot';  
defaultHeaders        = '';  
defaultNPlotsRow      = 6;
defaultNRowsPage      = 100;
saveFig               = false;
% 
p = inputParser;
addOptional(p,'cellIDStr',defaultCellIDStr,@(x) isstring(x) || iscell(x));
addParameter(p,'cmap',defaultCMap,@ischar);
addParameter(p,'nsteps',defaultNSteps,@isscalar);
addParameter(p,'plotsize',defaultPlotSize);
addParameter(p,'plotsep',defaultPlotSep);
addParameter(p,'offsetbase',defaultOffsetBase);
addParameter(p,'figname',defaultFigName,@ischar);
addParameter(p,'headers',defaultHeaders,@iscell);
addParameter(p,'nplots',defaultNPlotsRow,@isscalar);
addParameter(p,'nrows',defaultNRowsPage,@isscalar);
addParameter(p,'save',saveFig,@islogical);
parse(p,varargin{:});

% some sanity checcks should go here
% TO fix: GRID MODE IS BROKEN CURRENTLY - I THINK FIXED! 

%%
% gather plot tiling
if numel(data) == 1 
    noGridMode   = false;
else
    noGridMode   = true;
end
nPlots    = prod([size(data) length(data{1})]);
%
% open a waitbar
plotCount = 1;
hWait     = waitbar(0); 

nFigs = ceil(length(data{1})/p.Results.nrows);
rowCount = 1;
for n = 1:nFigs
    
    maps2plot = min([p.Results.nrows,length(data{1}) - (rowCount-1)]);
    
    % gather plot tiling
    if ~noGridMode
        nCols = min([maps2plot,p.Results.nplots]); % in case just 1 trial, we want to make a compact plot...
        nRows = ceil(maps2plot/nCols);
        %nRows = maps2plot;
    else
        nCols = size(data,2) * size(data,3); % ...if it's several trials, we plot nCells across trials
        nRows = maps2plot;
    end
    
    [axArray, hScroll] = scanpix.plot.multPlot([nRows nCols],'offsetbase',p.Results.offsetbase,'plotsep',p.Results.plotsep,'plotsize',p.Results.plotsize,'figname',p.Results.figname);
    hScroll.hFig.Visible = 'off';
    axCol = 1;
    axRow = 1;
    %
    
    for i = 1:maps2plot%nRows
        
        for k = 1:size(data,2)
            
            for j = 1:size(data,3)
                
                waitbar(plotCount/nPlots,hWait,'Making your precious figure, just bare with me!');
                
                if ~isempty(data{1,k,j}{i})
                    
                    
                    plotPeakRateFlag = false;
                    % plot
                    hAx = axArray{axRow,axCol};
                    
                    if strcmpi(type{j},'rate') || strcmpi(type{j},'spike') || strcmpi(type{j},'pos')
                        scanpix.plot.plotRateMap(data{1,k,j}{i},hAx,'colmap',p.Results.cmap,'nsteps',p.Results.nsteps);
                        % scale all axes to the largest possible map here? makes
                        % everything appear in proportion
                        if i == 1
                            axLims = max(cell2mat(cellfun(@size,vertcat(data{1,:,j}),'UniformOutput',false)),2);
                        end
                        set(hAx,'ydir','reverse','xlim',[0 axLims(2)],'ylim',[0 axLims(1)]);
                        plotPeakRateFlag = true;
                    elseif  strcmp(type{j},'dir')
                        scanpix.plot.plotDirMap(data{1,k,j}{i},hAx);
                        plotPeakRateFlag = true;
                    elseif  strcmp(type{j},'lin')
                        scanpix.plot.plotLinMaps( data{1,k,j}{i}, hAx);
                    elseif  strcmpi(type{j},'sacs')
                        plotPeakRateFlag = true;
                        mapSz = max(cell2mat(cellfun(@(x) size(x),data{1,k,j},'uni',0)));
                        imagesc(hAx,'CData',data{1,k,j}{i}); colormap(hAx,jet);
                        %                 axis(hAx,'square');
                        
                        set(hAx,'xlim',[0 mapSz(1)],'ylim',[0 mapSz(1)]);
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
                    
                    if i == nRows && noGridMode
                        if isempty(p.Results.headers)
                            title(hAx,['Trial_' num2str(k) '-' type{j}],'Interpreter','none');
                        else
                            title(hAx,[p.Results.headers{k,j}],'Interpreter','none');
                        end
                    end
                end
                
                axCol = axCol + 1;
                % bump indices
                if axCol > nCols
                    axRow = axRow + 1;
                    axCol = 1;
                    %rowCount = rowCount + 1;
                end

                plotCount = plotCount + 1;
            end
        end
    end
    %
    if p.Results.save
        scanpix.helpers.saveFigAsPDF([p.Results.figname '_' num2str(n)], [cd filesep]);
        close(hScroll.hFig);
    end
    %
    if exist('hScroll','var') && ishandle(hScroll.hFig)
        hScroll.hFig.Visible = 'on';
    end
end

close(hWait);


end

