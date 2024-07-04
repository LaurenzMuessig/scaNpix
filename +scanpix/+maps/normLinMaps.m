function [linMapsNormed, normInd] = normLinMaps(linMaps)
% normLinMaps - normalise a set of linear rate maps and sort by their peak
% position on the track
%
% package: scanpix.maps
%
% Syntax:
%       scanpix.maps.normLinMaps(linMaps)
%
% Inputs:
%    linMapsNormed - array of stacked linear rate maps  
% Outputs:
%    linMapsNormed - array of stacked normalised maps
%    normInd       - index for sorting maps
%
% see also:
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxMaps       = max(linMaps, [], 2, 'omitnan');
[~, maxInd]   = max( bsxfun(@eq, linMaps, maxMaps),[], 2, 'omitnan');
[~, normInd]  = sort( maxInd );
linMapsNormed = linMaps(normInd,:) ./ maxMaps(normInd); % ordered maps

end

