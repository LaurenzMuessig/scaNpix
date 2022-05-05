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
defaultPlotSep        = [50 30];  % pixel
defaultOffsetBase     = [60 50];  % pixel
defaultFigName        = 'multPlot';  
defaultHeaders        = '';  
defaultNPlotsRow      = 6;
% 
p = inputParser;
addOptional(p,'cellIDStr',defaultCellIDStr,@isstring);
addParameter(p,'cmap',defaultCMap,@ischar);
addParameter(p,'nsteps',defaultNSteps,@isscalar);
addParameter(p,'plotsize',defaultPlotSize);
addParameter(p,'plotsep',defaultPlotSep);
addParameter(p,'offsetbase',defaultOffsetBase);
addParameter(p,'figname',defaultFigName,@ischar);
addParameter(p,'headers',defaultHeaders,@iscell);
addParameter(p,'nplots',defaultNPlotsRow,@isscalar);
parse(p,varargin{:});

% some sanity checcks should go here

%%
% gather plot tiling
if numel(data) == 1 
%     nPlotsPerRow = p.Results.plotsize; % in case just 1 trial, we want to make a compact plot...
    noGridMode   = false;
    axArraySz    = [ceil(length(data{1})/p.Results.nplots) p.Results.nplots];
else
%     nPlotsPerRow = numel(data); % ...if it's several trials/plot types, we plot nCells across those
    noGridMode   = true;
    axArraySz    = [length(data{1})*size(data,1) numel(data)];
end
nPlots    = prod([size(data) length(data{1})]);
%
[axArray, hScroll] = scanpix.plot.multPlot(axArraySz,'plotsize',p.Results.plotsize,'plotsep',p.Results.plotsep,'figname',p.Results.figname);
% hScroll.hFig.Visible = 'off';

% open a waitbar
plotCount = 1;
hWait     = waitbar(0); 

for i = 1:length(data{1})
    
    for k = 1:size(data,2)
        
        for j = 1:size(data,3)
                      
            waitbar(plotCount/nPlots,hWait,'Making your precious figure, just bare with me!');
            
%             if isempty(data{1,k,j}{i})
%                 continue
%             end
            plotPeakRateFlag = false;
            % plot
            [b, a] = ind2sub(fliplr(axArraySz),plotCount);
            hAx = axArray{a,b};
            
            if strcmpi(type{j},'rate') || strcmpi(type{j},'spike') || strcmpi(type{j},'pos')
                scanpix.plot.plotRateMap(data{1,k,j}{i},hAx,'colmap',p.Results.cmap,'nsteps',p.Results.nsteps)
                plotPeakRateFlag = true;
            elseif  strcmp(type{j},'dir')
                scanpix.plot.plotDirMap(data{1,k,j}{i},hAx);
                plotPeakRateFlag = true;
            elseif  strcmp(type{j},'lin')
                 scanpix.plot.plotLinMaps( data{1,k,j}{i}, hAx);
            elseif  strcmpi(type{j},'sacs')
                plotPeakRateFlag = true;
                mapSz = size(data{1,k,j}{i});
                imagesc(hAx,'CData',data{1,k,j}{i}); colormap(hAx,jet);
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
                set(t,'Units','pixels','position',[-30 hAx.Position(4)/2],'String',p.Results.cellIDStr{i},'FontSize',8,'Interpreter','none' ); % harcoded text pos
            end
            
            if i == length(data{1}) && noGridMode
                if isempty(p.Results.headers)
                    title(hAx,['Trial_' num2str(k) '-' type{j}],'Interpreter','none');
                else
                    title(hAx,[p.Results.headers{k,j}],'Interpreter','none');
                end
            end
            plotCount = plotCount + 1;
        end
    end
end

close(hWait);

% hScroll.hFig.Visible = 'on';

end

