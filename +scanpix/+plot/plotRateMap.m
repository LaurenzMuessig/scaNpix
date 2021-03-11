function plotRateMap(map,varargin)
%
%
% parse input
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

% plot
[rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(map, p.Results.colmap, p.Results.nsteps); % bin rate map
% plot heat map
imagesc(hAx, rMapBinned);
caxis(hAx, [0 nanmax(rMapBinned(:))]);
colormap(hAx, cMapBinned);
axis(hAx,'off');

end

