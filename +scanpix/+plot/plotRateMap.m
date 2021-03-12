function plotRateMap(map,varargin)
% plotRateMap - plot a standard rate map
% package: scanpix.plot
%
%  Usage:   scanpix.plot.plotRateMap( dirMap ) 
%           scanpix.plot.plotRateMap( dirMap, hAx )
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
defaultColMap = 'jet';
defaultNSteps = 11;
defaulthAx    = [];

p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'colmap',defaultColMap,@ischar);
addParameter(p,'nsteps',defaultNSteps,@isscalar);
parse(p,varargin{:});

if isempty(p.Results.ax)
    hAx = axes;
else
    hAx = p.Results.ax;
end

%% plot
[rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(map, p.Results.colmap, p.Results.nsteps); % bin rate map
% plot heat map
imagesc(hAx, rMapBinned);
caxis(hAx, [0 nanmax(rMapBinned(:))]);
colormap(hAx, cMapBinned);
axis(hAx,'off');

end

