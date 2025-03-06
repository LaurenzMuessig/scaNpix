function plotRateMap(map,ax,options)
% plotRateMap - plot a standard rate map
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotRateMap( rate map ) 
%           scanpix.plot.plotRateMap( rate map, hAx )
%           scanpix.plot.plotRateMap( __ ,'name',value,.... )
%
%  Inputs:  
%           map      - rate map
%           varargin - optional inputs
%                    - axis handle and/or
%                    - Name-Value pairs for 'colmap' (Matlab color map)
%                      and/or 'nsteps' (how many steps in colour map)
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parse input
arguments
    map {mustBeNumeric}
    ax  {ishghandle(ax, 'axes')} = axes;
    options.colmap = 'jet' ;
    options.nsteps (1,1) {mustBeNumeric} = 11;
    options.interp (1,1) {mustBeNumericOrLogical} = false;
    options.cmapEdge (1,2) {mustBeNumeric} = [0 max(map(:),[],'omitnan')];
end

%%
if isempty(map);  return; end

%% plot
% bin rate map
[rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(map, 'colmap',options.colmap, 'nsteps', options.nsteps, 'cmapEdge', options.cmapEdge); % bin rate map
%
if options.interp
    rMapBinned = scanpix.maps.interpMap(rMapBinned); 
end

% plot heat map
imagesc(ax,'CData',rMapBinned,[0 size(cMapBinned,1)]);
colormap(ax, cMapBinned);
axis(ax,'off');
%
mapSz = size(map);
set(ax,'ydir','reverse','xlim',[0.5 mapSz(2)+0.5],'ylim',[0.5 mapSz(1)+0.5]);
% set(ax,'xlim',[0.5 mapSz(2)+0.5],'ylim',[0.5 mapSz(1)+0.5]);
xLabelStruct = get(ax,'xlabel');
xLabelStruct.VerticalAlignment = 'bottom'; 
xLabelStruct.String = ['peakFR=' num2str(max(map(:),[],'omitnan'),'%.1f')];
xLabelStruct.Visible = 'On';
%
if mapSz(1) == mapSz(2)
    axis(ax,'square');
end

end

