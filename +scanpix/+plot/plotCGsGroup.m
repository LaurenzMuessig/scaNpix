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
    if isempty(p.Results.maps); error('scaNpix::plot::plotCGsGroup: You need to supply maps if you want to plot them. Not sure why I have to tell you that...'); end
end

if length(spikeTimes) < 2
    warning('scaNpix::plot::plotCGsGroup: You need to supply spike times for at least 2 cells to make this figure. It''s a comparison figure, duh.');
    return
end

%% plot
%wait bar
hWait         = waitbar(0); 
plotCount     = 1;

[axArr, hScroll] = scanpix.plot.multPlot([length(spikeTimes) length(spikeTimes)],'plotsize',p.Results.plotsize,'plotsep',p.Results.plotsep,'offsetbase',p.Results.offsetbase,'figname',p.Results.figname);
nPlots           = numel(axArr)/2;

hScroll.hFig.Visible = 'off';

for i = 1:length(spikeTimes)
    % wait bar
    waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
    % 
    if p.Results.plotmaps
        scanpix.plot.plotRateMap(p.Results.maps{i},axArr{i,i});
        % text(axArr{i,i},-12,0.45*max(get(axArr{i,i},'ylim')),p.Results.cellIDStr{i},'Interpreter','none');
    else
        % plot AC
        scanpix.analysis.spk_crosscorr(spikeTimes{i},'AC',binSz,lag,trialDur,'plot',axArr{i,i}); % autocorr
        if i ~= 1
            set(axArr{i,i},'xtick',[-lag 0 lag],'xticklabel',{});
        else
            set(axArr{i,i},'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
        end

    end
    % yAxlim    = get(axArr{i,i},'ylim');
    % plot headers
    text(axArr{i,i},-0.65,0.5,p.Results.cellIDStr{i},'Units','normalized','Interpreter','none');
    text(axArr{i,i},0.2,1.1,p.Results.cellIDStr{i},'Units','normalized','Interpreter','none');
    %
    plotCount = plotCount + 1;
    for j = i+1:length(spikeTimes)

        waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
        % plot cross corr
        scanpix.analysis.spk_crosscorr(spikeTimes{i},spikeTimes{j},binSz,lag,trialDur,'plot',axArr{i,j}); % crosscorr
        yAxlim = get(axArr{i,j},'ylim');
         % plot comparison ID
        % text(axArr{i,j},-lag-0.0025,1.15*max(yAxlim),[p.Results.cellIDStr{i} ' v ' p.Results.cellIDStr{j}],'Interpreter','none','color','r');
        if i~=1
            set(axArr{i,j},'xtick',[-lag 0 lag],'xticklabel',{''});
        else
            set(axArr{i,j},'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
        end
        if j == length(spikeTimes)
            text(axArr{i,j},1.1*lag,0.5*max(yAxlim),p.Results.cellIDStr{i},'Interpreter','none');
        end
        plotCount  = plotCount + 1;
    end
end

% if p.Results.plotmaps
%     scanpix.plot.plotRateMap(p.Results.maps{i},axArr{end,end});
%     text(axArr{end,end},max(get(axArr{end,end},'xlim'))+3,0.45*max(get(axArr{end,end},'ylim')),p.Results.cellIDStr{i},'Interpreter','none');
% else
    text(axArr{i,j},1.1*lag,0.5*max(yAxlim),p.Results.cellIDStr{i},'Interpreter','none');
% end
% remove empty axes
remInd = cellfun(@(x) isempty(get(x, 'Children')),axArr);
delete([axArr{remInd}])
%
close(hWait);
hScroll.hFig.Visible = 'on';
end


