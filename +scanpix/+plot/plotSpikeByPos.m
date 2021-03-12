function plotSpikeByPos(xy,spikeT,varargin)
% plotSpikeByPos - plot location of spikes of a given cell on top of path
% of animal during trial
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotSpikeByPos( xy, spikeT) 
%           scanpix.plot.plotSpikeByPos( xy, spikeT, hAx )
%           scanpix.plot.plotSpikeByPos( __ ,'srate',value,)
%
%  Inputs:  
%           xy      - xy positions (nSamples-by-2)
%           spikeT  - spike times (in seconds)
%           varargin - optional inputs
%                    - axis handle and/or
%                    - Name-Value pairs for 'srate' (positional sampling rate)
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% parse input
defaultPosSR = 50;
defaulthAx    = [];
% 
p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'srate',defaultPosSR,@isscalar);
parse(p,varargin{:});
% 
if isempty(p.Results.ax)
    hAx = axes;
else
    hAx = p.Results.ax;
end

%% plot
spikeInd = ceil(spikeT .* p.Results.srate); % bin spikjes

% plot
plot(hAx,xy(:,1),xy(:,2),'k-','linewidth',0.5);
hold(hAx, 'on');
plot(hAx,xy(spikeInd,1),xy(spikeInd,2),'rs','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',2.5);
hold(hAx, 'off');
axis(hAx, 'off');

end

