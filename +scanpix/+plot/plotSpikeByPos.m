function plotSpikeByPos(xy,spikeT,ax, options)
% plotSpikeByPos - plot location of spikes of a given cell on top of path
% of animal during trial
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotSpikeByPos( xy, spikeT) 
%           scanpix.plot.plotSpikeByPos( xy, spikeT, hAx )
%           scanpix.plot.plotSpikeByPos( __ ,'srate',value,)
%
%  Inputs:  
%           xy       - xy positions (nSamples-by-2)
%           spikeT   - spike times (in seconds)
%           varargin - optional inputs - see parse inputs section
%                    - axis handle and/or
%                    - Name-Value pairs for 'srate' (positional sampling
%                      rate), 'msize' (marker size) and/or 'axvis' (flag for
%                      axis on/off)
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%% To-Do:
% Add using sample times for binning spikes (or make binned spike times the
% input?

%% parse input
arguments
    xy {mustBeNumeric}
    spikeT {mustBeNumeric}
    ax  {ishghandle(ax, 'axes')} = axes;
    options.srate (1,1) {mustBeNumeric} = 50 ;
    options.msize (1,1) {mustBeNumeric} = 2.5;
    options.mcol   = 'r';
    options.axvis (1,1) {mustBeNumericOrLogical} = false;
end

tf = ishold(ax);

%% plot
spikeInd = ceil(spikeT .* options.srate); % bin spikes

% plot
% plot(ax,xy(:,1),xy(:,2),'k-',xy(spikeInd,1),xy(spikeInd,2),'s');
% if isempty(spikeInd)
%     ax.Children(1).LineWidth = 0.5;
% else
%     ax.Children(1).LineStyle = 'none'; ax.Children(1).MarkerFaceColor = options.mcol; ax.Children(1).MarkerEdgeColor = 'none'; ax.Children(1).MarkerSize = options.msize; % spike pos  
%     ax.Children(2).LineWidth = 0.5; % path
% end

plot(ax,xy(:,1),xy(:,2),'k-','linewidth',0.5);
if ~isempty(spikeInd)
    if ~tf; hold(ax, 'on'); end
    plot(ax,xy(spikeInd,1),xy(spikeInd,2),'s','LineStyle','none','MarkerFaceColor',options.mcol,'MarkerEdgeColor','none','MarkerSize',options.msize);
    if ~tf; hold(ax, 'off'); end
end

if ~options.axvis; axis(ax, 'off'); end
end

