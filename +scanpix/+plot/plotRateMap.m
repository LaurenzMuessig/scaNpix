function plotRateMap(map,varargin)
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
defaultColMap   = 'jet';
defaultNSteps   = 11;
defaulthAx      = [];
defaultCMapEdge = [];

p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'colmap',defaultColMap,@ischar);
addParameter(p,'nsteps',defaultNSteps,@isscalar);
addParameter(p,'cmapEdge',defaultCMapEdge);
parse(p,varargin{:});

if isempty(p.Results.ax)
    hAx = axes;
else
    hAx = p.Results.ax;
end

% if strcmpi(p.Results.colmap,'poulter')
%     nSteps = 8;  % Steve's map has 8 fixed steps, so should ignore any dynamic setting here
% else
%     nSteps = p.Results.nsteps;
% end

%% plot
if isempty(p.Results.cmapEdge)
    [rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(map, 'colmap',p.Results.colmap, 'nsteps', p.Results.nsteps); % bin rate map
else
    [rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(map, 'colmap',p.Results.colmap, 'nsteps', p.Results.nsteps, 'cmapEdge', p.Results.cmapEdge); % bin rate map
end
% plot heat map
imagesc(hAx,'CData',rMapBinned,[0 size(cMapBinned,1)]);

colormap(hAx, cMapBinned);
axis(hAx,'off');

% imagesc(hAx, rMapBinned);
% caxis(hAx, [0 p.Results.nsteps+1]);
% colormap(hAx, cMapBinned);

end

