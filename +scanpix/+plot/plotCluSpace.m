function plotCluSpace(waveForms,scaleMax,options)
% Callback function for 'plotCLUSpace' button 
% This function will plot/reconstruct the cluster space (as in Tint) 
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TO-Do:
% add callback so different clusters can be highlighted by user
%% parse input
arguments
    waveForms {mustBeA(waveForms,'cell')} 
    scaleMax {mustBeNumeric} 
    options.figName (1,:) {mustBeText} = 'scaNpix::CluSpace';
    options.cellLabels {mustBeA(options.cellLabels,'cell')}  = strcat({'cell_'},num2str((1:length(waveForms))'));
end

%%

% peak-trough amplitudes
funH          = @(a) squeeze(max(a,[],2) - min(a,[],2)); % handle fo cellfun call (peak to trough ampl)
amplitudes    = cellfun(funH,waveForms,'UniformOutput', 0);
channelComp   = nchoosek(1:4,2); % could also just hard code
scaleMax      = scaleMax(channelComp);

% make figure
nRows         = ceil(mod(size(channelComp,1),4));
nCols         = ceil(size(channelComp,1) / nRows);
axArray       = scanpix.plot.multPlot([nRows,nCols],'figname',options.figName,'plotsize',[250 250]);
axArray       = axArray(:);
% 
col           = scanpix.fxchange.cbrewer('qual', 'Set1', max([length(amplitudes), 3]), 'PCHIP' ); % cbrewer makes nice colormaps
for i = 1:size(channelComp,1)
    
    hold(axArray{i},'on');
    % plot cluster space
    for j = 1:length(amplitudes)
        if ~isempty(amplitudes{j})
            scatter(axArray{i},amplitudes{j}(:,channelComp(i,2)),amplitudes{j}(:,channelComp(i,1)),6,'filled','markerfacecolor',col(j,:),'markeredgecolor','none','markerfacealpha',0.6);
        end
    end
    % loop again to plot cluster IDs (otherwise they might not be visible in a busy space)
    for j = 1:length(amplitudes)
        %%% DEAL WITH TEXT COLOUR %%%%
        if ~isempty(amplitudes{j})
            text(axArray{i},mean(amplitudes{j}(:,channelComp(i,2))),mean(amplitudes{j}(:,channelComp(i,1))),options.cellLabels{j},'color',[0 0 0],'FontSize',12,'FontWeight','bold','HorizontalAlignment','center'); % p
        end
        % plot at centroid
    end
    hold(axArray{i},'off');
    % some more axis formatting
    set(axArray{i},'xlim',[0 2*scaleMax(channelComp(i,2))],'ylim',[0 2*scaleMax(channelComp(i,1))],'xtick',[],'xticklabel',{''},'ytick',[],'yticklabel',{''},'box','on');
    axis(axArray{i},'square');
    xlabel(axArray{i},['channel_' num2str(channelComp(i,2))],'interpreter','none','fontsize',12,'fontweight','bold'); 
    ylabel(axArray{i},['channel_' num2str(channelComp(i,1))],'interpreter','none','fontsize',12,'fontweight','bold');
end

end

