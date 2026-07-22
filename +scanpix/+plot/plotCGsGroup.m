function plotCGsGroup(spikeTimes,binSz,lag,trialDur,options)
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


%%
arguments
    spikeTimes {mustBeA(spikeTimes,'cell')}
    binSz (1,1) {mustBeNumeric}
    lag (1,1) {mustBeNumeric}
    trialDur (1,1) {mustBeNumeric}
    options.maps {mustBeA(options.maps,'cell')} = {};
    options.cellIDStr {mustBeA(options.cellIDStr,'cell')} = strcat({'cell_'},num2str((1:length(spikeTimes))'));
    options.figname (1,:) {mustBeText} = 'scaNpix::CG_overview';
    options.plotsize (1,2) {mustBeNumeric} = [75 75];
    options.plotsep (1,2) {mustBeNumeric} = [15 20];
    options.offset (1,2) {mustBeNumeric} = [50 40];
    options.plotmaps (1,1) {mustBeNumericOrLogical} = false;
    % options.groupInd (1,1) {mustBeNumeric}
end

%%
if length(spikeTimes) < 2
    warning('scaNpix::plot::plotCGsGroup: You need to supply spike times for at least 2 cells to make this figure. It''s a comparison figure, duh.');
    return
end

if options.plotmaps && isempty(options.maps)
    warning('scaNpix::plot::plotCGsGroup: You need to supply maps if you want to plot them. Not sure why I have to tell you that...');
    options.plotmaps = false;
end



%% plot
%wait bar
hWait         = waitbar(0); 
plotCount     = 1;

[axArr, hScroll] = scanpix.plot.multPlot([length(spikeTimes) length(spikeTimes)],'plotsize',options.plotsize,'plotsep',options.plotsep,'offset',options.offset,'figname',options.figname);
nPlots           = numel(axArr)/2;

hScroll.hFig.Visible = 'off';

for i = 1:length(spikeTimes)
    % wait bar
    waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
    % 
    if options.plotmaps
        scanpix.plot.plotRateMap(options.maps{i},axArr{i,i});
        % text(axArr{i,i},-12,0.45*max(get(axArr{i,i},'ylim')),options..cellIDStr{i},'Interpreter','none');
    else
        % plot AC
        scanpix.analysis.spk_crosscorr(spikeTimes{i},'AC',binSz,lag,trialDur,'plot',axArr{i,i}); % autocorr
        if i ~= 1
            set(axArr{i,i},'xtick',[-lag 0 lag],'xticklabel',{});
        else
            set(axArr{i,i},'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
        end

    end
    % yAxlim = get(axArr{i,i},'ylim');
    % plot headers
    text(axArr{i,i},-0.65,0.5,options.cellIDStr{i},'Units','normalized','Interpreter','none');
    text(axArr{i,i},0.2,1.1,options.cellIDStr{i},'Units','normalized','Interpreter','none');
    %
    plotCount = plotCount + 1;
    for j = i+1:length(spikeTimes)

        waitbar(plotCount/nPlots,hWait,'Plotting Correlograms, just bare with me!');
        % plot cross corr
        scanpix.analysis.spk_crosscorr(spikeTimes{i},spikeTimes{j},binSz,lag,trialDur,'plot',axArr{i,j}); % crosscorr
        % plot comparison ID
        % text(axArr{i,j},-lag-0.0025,1.15*max(yAxlim),[options..cellIDStr{i} ' v ' options..cellIDStr{j}],'Interpreter','none','color','r');
        if i~=1
            set(axArr{i,j},'xtick',[-lag 0 lag],'xticklabel',{''});
        else
            set(axArr{i,j},'xtick',[-lag 0 lag],'xticklabel',[-lag*1000,0,lag*1000]);
        end
        %
        yAxlim = get(axArr{i,j},'ylim');
        if j == length(spikeTimes)
            text(axArr{i,j},1.1*lag,0.5*max(yAxlim),options.cellIDStr{i},'Interpreter','none');
        end
        plotCount  = plotCount + 1;
    end
end
%
text(axArr{i,j},1.1*lag,0.5*max(yAxlim),options.cellIDStr{i},'Interpreter','none');
% remove empty axes
scanpix.plot.cleanAxMultiPlot(axArr,'del');
%
close(hWait);
hScroll.hFig.Visible = 'on';
end

% if options..plotmaps
%     scanpix.plot.plotRateMap(options..maps{i},axArr{end,end});
%     text(axArr{end,end},max(get(axArr{end,end},'xlim'))+3,0.45*max(get(axArr{end,end},'ylim')),options..cellIDStr{i},'Interpreter','none');
% else
    
% end

