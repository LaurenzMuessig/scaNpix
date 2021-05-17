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
[x,y] = pol2cart(thetas, rhos);
plot(hAx, [x; x(1)], [y; y(1)],'k','linewidth',2);
hold(hAx,'on');
% plot mean vector and format axis
[x,y] = pol2cart(meanDir,meanR);
plot(hAx,[0 x],[0 y],'r-','linewidth',2); % mean vector
plot(hAx,[-1 1; 0 0]',[0 0; -1 1]','k-','linewidth',1.5);
hold(hAx,'off');
axis(hAx,'square');
end

