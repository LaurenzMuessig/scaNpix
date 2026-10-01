function [meanR,meanDir,thetas,rhos] = rayleighVect(dirMap)
% rayleighVect - calculate Rayleigh vector length and direction from a 
% directional rate map. Also ouput angles and radii from underlying data as
% convenient in case you want to also plot the data 
% package: scanpix.analysis
%
%  Usage:   scanpix.analysis.rayleighVect( dirMap )
%
%  Inputs:  
%           dirMap - directional rate map
%
% Outputs: 
%           meanR   - rayleigh vector length
%           meanDir - rayleigh vector direction (rad, 0 <= meanDir < 2*pi; NaN for empty maps)
%           thetas  - binned angles
%           rhos    - magnitude for each bin (norm. firing rate)
%
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    dirMap (:,1) {mustBeNumeric}

end
%%
binSz   = 2*pi / length(dirMap); % binSz

thetas  = linspace(binSz/2, 2*pi-binSz/2, length(dirMap))'; % binned angles
rhos    = dirMap ./ max(dirMap(:),[],'omitnan');            % normalised rates

% get mean vector
rVect   = sum(rhos .* exp(1i*thetas), 'omitnan');           % rayleigh vector, complex (rate weighted sum of unit vectors)
meanR   = abs(rVect)/sum(rhos,'omitnan');                   % length of resultant vector normalised by sum of weights
meanDir = mod(angle(rVect), 2*pi);                          % direction of resultant vector (in rad, 0 <= meanDir < 2*pi)
if meanDir == 2*pi; meanDir = 0; end                        % mod can round tiny negative angles up to exactly 2*pi

% no spikes / empty map - direction undefined
if ~isfinite(meanR)
    meanDir = NaN;
end

end

