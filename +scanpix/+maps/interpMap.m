function map_interp = interpMap(map, nBins)
% interpMap - interpolate a 2D map, e.g. rate map or spatial AC
% package: scanpix.map
%
% Interpolate a rate map to make the coloring appear smoother
%
%  Usage:   rate_map_interp = scanpix.plot.interpMap( rate_map ) 
%
%  Inputs:  
%           rate_map        - rate map (any 2D map type really) 
%
%  Outputs:
%           rate_map_interp - interpolated rate map
%
% LM 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    map {mustBeNumeric}
    nBins (1,1) {mustBeNumeric} = 1000;
end

%%
bins4InterpX = linspace(1,size(map,2),nBins);
bins4InterpY = linspace(1,size(map,1),nBins);

%% interpolate NaNs in rate map as average from surrounding non-NaN bins (so NaNs don't spread during interpolations)
nanInd            = isnan(map);
mapNoNans         = map;
map(nanInd)       = 0;
visMask           = ones(size(map));
visMask(nanInd)   = 0;
visMask_sm        = imfilter(visMask, ones(3));
tmpMap            = imfilter(map, ones(3)) ./ visMask_sm;
mapNoNans(nanInd) = tmpMap(nanInd);

%% interpolate
[X, Y]            = meshgrid(1:size(map,2), 1:size(map,1));
[X2, Y2]          = meshgrid(bins4InterpX, bins4InterpY);
map_interp        = interp2(X, Y, mapNoNans, X2, Y2, 'linear');

end

