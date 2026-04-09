function [vMapC, vMapExtentMask] = circulariseVMap( vMap, angBinInDeg, binInterpF )
% Convert to vMap to polar coordinates, and resample to an XY grid. Hence gives an 'image'-friendly
% array of circularised vector map.

%%
arguments
    vMap {mustBeNumeric} 
    angBinInDeg (1,1) {mustBeNumeric} = 6;
    binInterpF (1,1) {mustBeNumeric} = 5;  % Sets how many bins are in the *square image* of circ vMap (units=relative to distance bin of cartesian vMap). Creates oversampling at the edge, but allows more visible detail in the centre.
end

%%
if isempty(vMap);   vMapC = [];   return;   end

%
angBins                  = ((angBinInDeg:angBinInDeg:360)./360).*2.*pi;
dBins                    = 1:size(vMap,1);
[phiBinsArr,dBinsArr]    = meshgrid( angBins, dBins );
[polBinsInX, polBinsInY] = pol2cart( phiBinsArr, dBinsArr );
[resampX, resampY]       = meshgrid( -max(dBins):(1/binInterpF):max(dBins) );
% First resample the firing rates.
resampR                  = griddata( polBinsInX(:), polBinsInY(:), vMap(:), resampX(:), resampY(:) );
vMapC                    = reshape( resampR, size(resampX) );
% Also resample the vect map itself (all set to 0), so we can mark unvisited vMap differently to circ vMap surround.
vMapDum                  = zeros(size(vMap));
resampDum                = griddata( polBinsInX(:), polBinsInY(:), vMapDum(:), resampX(:), resampY(:) );
vMapExtentMask           = reshape( resampDum, size(resampX) );

end