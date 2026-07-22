function plotSpeedMap(map,options)
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

%%
arguments
    map {mustBeNumeric}
    options.ax  {ishghandle(options.ax, 'axes')} = axes;
    options.maxspeed (1,1) {mustBeNumeric} = 40;
end

%% plot
confInt = [map(:,2) - map(:,3), map(:,2) + map(:,3); NaN NaN];
xVals   = 1:size(map(:,2),1)+1;
plot(options.ax,xVals(1:end-1)',map(:,2),'r-',[xVals'; xVals'],confInt(:),'r:');
%
options.ax.Children(1).LineWidth = 1;
options.ax.Children(2).LineWidth = 2;

% format axis
if sum(map(:,2)) ~= 0
    set(options.ax,'ylim',[min([0;confInt(:)]) max(confInt(:))],'ytick',[0 max(map(:,2),[],'omitnan')],'YTickLabel',{'0' sprintf('%.1f',max(map(:,2),[],'omitnan'))},'xlim',[0 length(map(:,1))+1],'xtick',[0 length(map(:,1))+1],'xticklabel',[0 round(options.maxspeed)]);
end
%
ylabel(options.ax,'Firing Rate (Hz)');
xlabel(options.ax,'Running Speed (cm/s)','VerticalAlignment','middle');
axis(options.ax,'square');

end

