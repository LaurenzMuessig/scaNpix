function plotDirMap(dirMap,hAx)
% plotDirMap - plot a directional rate map as a standard polar plot
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotDirMap( dirMap ) 
%           scanpix.plot.plotDirMap( dirMap, hAx )
%
%  Inputs:  
%           dirMap - directional rate map
%           hAx    - axis handle (optional)
%
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    hAx = axes;
end

%
[meanR,meanDir,thetas,rhos] = scanpix.analysis.rayleighVect(dirMap);

% plot as line rather than 'polarplot' as I don't like the circle outline
% of the polarplots in Matlab
plot(hAx,[-1 1; 0 0]',[0 0; -1 1]','k-','linewidth',1.5); % axis
hold(hAx,'on');
[x,y] = pol2cart(thetas, rhos);
plot(hAx, [x; x(1)], [y; y(1)],'k','linewidth',2);
% plot mean vector 
[x,y] = pol2cart(meanDir,meanR);
plot(hAx,[0 x],[0 y],'r-','linewidth',2);
% format axis
hold(hAx,'off');
axis(hAx,'square','off');
% axis(hAx,'off');
end

