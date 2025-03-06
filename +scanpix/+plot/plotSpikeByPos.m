function plotSpikeByPos(xy,spikeT,ax,options)
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
    options.msize (1,1) {mustBeNumeric} = 10;
    options.mcol   = 'r';
    options.axvis (1,1) {mustBeNumericOrLogical} = false;
end

tf = ishold(ax);

%% plot
spikeInd = ceil(spikeT .* options.srate); % bin spikes

% plot
plot(ax,xy(:,1),xy(:,2),'k-','linewidth',0.5);
if ~isempty(spikeInd)
    if ~tf; hold(ax, 'on'); end
    % plot(ax,xy(spikeInd,1),xy(spikeInd,2),'s','LineStyle','none','MarkerFaceColor',options.mcol,'MarkerEdgeColor','none','MarkerSize',options.msize);
    scatter(ax,xy(spikeInd,1),xy(spikeInd,2),options.msize,options.mcol,'Filled'); %'MarkerEdgeColor','none');

    if ~tf; hold(ax, 'off'); end
end
set(ax,'ydir','reverse');

if ~options.axvis; axis(ax, 'off'); end
end

