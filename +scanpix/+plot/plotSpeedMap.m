function plotSpeedMap(map,varargin)
% plotSpeedMap - plot speed map
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotRateMap( rate map ) 
%           scanpix.plot.plotRateMap( rate map, hAx )
%           scanpix.plot.plotRateMap( __ ,'name',value,.... )
%
%  Inputs:  
%           map      - speed map
%           varargin - optional inputs
%                    - axis handle and/or
%                    - Name-Value pair for 'maxspeed'
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parse input
defaulthAx      = [];
defaultMaxSpeed = [];

p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'maxspeed',defaultMaxSpeed,@isscalar);
parse(p,varargin{:});

if isempty(p.Results.ax)
    hAx = axes;
else
    hAx = p.Results.ax;
end

if isempty(p.Results.maxspeed)
    maxSpeed = map(end,1);
else
    maxSpeed = p.Results.maxspeed;
end


%% plot
confInt = [map(:,2) + map(:,3:4); NaN NaN];
xVals = 1:size(map(:,2),1)+1;
plot(hAx,xVals(1:end-1)',map(:,2),'r-',[xVals'; xVals'],confInt(:),'r:');
hAx.Children(1).LineWidth = 1;
hAx.Children(2).LineWidth = 2;

% format axis
if sum(map(:,2)) ~= 0
    set(hAx,'ylim',[min([0;confInt(:)]) max(confInt(:))],'ytick',[0 nanmax(map(:,2))],'YTickLabel',{'0' sprintf('%.1f',nanmax(map(:,2)))},'xlim',[0 length(map(:,1))+1],'xtick',[0 length(map(:,1))+1],'xticklabel',[0 round(maxSpeed)]);
end

ylabel(hAx,'Firing Rate (Hz)');
xlabel(hAx,'Running Speed (cm/s)','VerticalAlignment','middle');
axis(hAx,'square');



end

