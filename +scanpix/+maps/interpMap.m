function map_interp = interpMap(map)
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
end

%%
[X, Y]                = meshgrid(1:size(map,2), 1:size(map,1));
[X2, Y2]              = meshgrid(1:0.01:size(map,2), 1:0.01:size(map,1));
nanInterp             = interp2(X, Y, double(isnan(map)), X2, Y2, 'linear') == 1;
map_interp            = interp2(X, Y, map, X2, Y2, 'linear');
map_interp(nanInterp) = NaN;

end

