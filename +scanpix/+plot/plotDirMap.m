function plotDirMap(dirMap,ax)
% plotDirMap - plot a directional rate map as a standard polar plot
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotDirMap( dirMap ) 
%           scanpix.plot.plotDirMap( dirMap, hAx )
%
%  Inputs:  
%           dirMap - directional rate map
%           ax     - axis handle (optional)
%
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
arguments
    dirMap {mustBeNumeric}
    ax     {ishghandle(ax, 'axes')} = axes;
end

%%
[meanR,meanDir,thetas,rhos] = scanpix.analysis.rayleighVect(dirMap);

%% plot
% plot as line rather than 'polarplot' as I don't like the circle outline
% of the polarplots in Matlab
plot(ax,[-1 1; 0 0]',[0 0; -1 1]','k-','linewidth',1.5); % axis
hold(ax,'on');
[x,y] = pol2cart(thetas, rhos);
plot(ax, [x; x(1)], [y; y(1)],'k','linewidth',2);
% plot mean vector 
[x,y] = pol2cart(meanDir,meanR);
plot(ax,[0 x],[0 y],'r-','linewidth',2);
% format axis
hold(ax,'off');
axis(ax,'square','off');
xLabelStruct = get(ax,'xlabel');
xLabelStruct.VerticalAlignment = 'bottom'; 
xLabelStruct.String = ['peakFR=' num2str(max(dirMap,[],'omitnan'),'%.1f')];
xLabelStruct.Visible = 'On';
end

