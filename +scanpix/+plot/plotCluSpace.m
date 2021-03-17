function plotCluSpace(waveForms,scaleMax,cell_IDs)
% Callback function for 'plotCLUSpace' button 
% This function will plot/reconstruct the cluster space (as in Tint) 
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TO-Do:
% add callback so different clusters can be highlighted by user

figSzPix = [550 500];

%% prepare
if nargin < 3
    cell_IDs = strcat({'cell_'},num2str((1:length(waveForms))'));
end

% peak-trough amplitudes
funH          = @(a) squeeze(max(a,[],2) - min(a,[],2)); % handle fo cellfun call (peak to trough ampl)
amplitudes    = cellfun(funH,waveForms,'UniformOutput', 0);
channelComp   = nchoosek(1:4,2); % could also just hard code
scaleMax      = scaleMax(channelComp);
% make figure
screenSz      = get(0,'screensize');
figure('units','pixel','position',[0.1*screenSz(3) 0.1*screenSz(4) figSzPix]);

offset        = [0 0];
col           = scanpix.fxchange.cbrewer('qual', 'Set1', max([length(amplitudes), 3]), 'PCHIP' ); % cbrewer makes nice colormaps
for i = 1:size(channelComp,1)
    
    axes('units','normalized','position',[0.05+offset(1) 0.6-offset(2) 0.275 0.25]);
    hold on
    % plot cluster space
    for j = 1:length(amplitudes)
        scatter(amplitudes{j}(:,channelComp(i,2)),amplitudes{j}(:,channelComp(i,1)),6,'filled','markerfacecolor',col(j,:),'markeredgecolor','none','markerfacealpha',0.6);
    end
    % loop again to plot cluster IDs (otherwise they might not be visible in a busy space)
    for j = 1:length(amplitudes)
        %%% DEAL WITH TEXT COLOUR %%%%
        text(mean(amplitudes{j}(:,channelComp(i,2))),mean(amplitudes{j}(:,channelComp(i,1))),cell_IDs{j},'color',[0 0 0],'FontSize',12,'FontWeight','bold','HorizontalAlignment','center'); % plot at centroid
    end
    hold off
    % some more axis formatting
    set(gca,'xlim',[0 2*scaleMax(channelComp(i,2))],'ylim',[0 2*scaleMax(channelComp(i,1))],'xtick',[],'xticklabel',{''},'ytick',[],'yticklabel',{''},'box','on');
    axis square
    xlabel(['channel_' num2str(channelComp(i,2))],'interpreter','none','fontsize',12,'fontweight','bold'); 
    ylabel(['channel_' num2str(channelComp(i,1))],'interpreter','none','fontsize',12,'fontweight','bold');
    if i ~= 3
        offset(1) = offset(1) + 0.33;
    else
        offset = [0 0.4];
    end
end

end

