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
%    linMapsNormed - cell array of linear rate maps  
% Outputs:
%    linMapsNormed - cell array of normalised maps
%    normInd       - index for sorting maps
%
% see also:
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linMapsNormed = cell(size(linMaps));

for i = 1:length(linMaps)
    maxMaps          = nanmax(linMaps{i}, [], 2);
    [~, maxInd]      = nanmax( bsxfun(@eq, linMaps{i}, maxMaps),[], 2 );
    [~, normInd]     = sort( maxInd );
    linMapsNormed{i} = linMaps{i}(normInd,:) ./ maxMaps(normInd); % ordered maps
end


end

