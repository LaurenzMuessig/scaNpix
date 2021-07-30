function plotSpikeByPos(xy,spikeT, varargin)
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
defaultPosSR    = 50;
defaulthAx      = [];
defaultMarkerSz = 2.5;   % marker size'
axisVisible     = false; % axis on/off
% 
p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'srate',defaultPosSR,@isscalar);
addParameter(p,'msize',defaultMarkerSz,@isscalar);
addParameter(p,'axvis',axisVisible,@islogical);
parse(p,varargin{:});
% 
if isempty(p.Results.ax)
    hAx = axes;
else
    hAx = p.Results.ax;
end

%% plot
spikeInd = ceil(spikeT .* p.Results.srate); % bin spikes

% plot
plot(hAx,xy(:,1),xy(:,2),'k-',xy(spikeInd,1),xy(spikeInd,2),'rs')
if isempty(spikeInd)
    hAx.Children(1).LineWidth = 0.5;
else
    hAx.Children(1).LineStyle = 'none'; hAx.Children(1).MarkerFaceColor = 'r'; hAx.Children(1).MarkerEdgeColor = 'none'; hAx.Children(1).MarkerSize = p.Results.msize; % spike pos  
    hAx.Children(2).LineWidth = 0.5; % path
end
% 'linewidth',0.5);
% plot(hAx,xy(:,1),xy(:,2),'k-','linewidth',0.5);
% hold(hAx, 'on');
% plot(hAx,xy(spikeInd,1),xy(spikeInd,2),'rs','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',p.Results.msize);
% hold(hAx, 'off');
if ~p.Results.axvis; axis(hAx, 'off'); end
end

