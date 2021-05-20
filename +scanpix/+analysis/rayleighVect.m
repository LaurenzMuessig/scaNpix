function [meanR,meanDir,thetas,rhos] = rayleighVect(dirMap)
% rayleighVect - calculate Rayleigh vector from a directional rate maps
% package: scanpix.analysis
%
%  Usage:   scanpix.analysis.rayleighVect( dirMap )
%
%  Inputs:  
%           dirMap - directional rate map
%
% Outputs: 
%           meanR   - rayleigh vector length
%           meanDir - rayleigh vector direction
%           thetas  - binned angles
%           rhos    - magnitude for each bin (norm. firing rate)
%
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
binSz  = 2*pi / length(dirMap); % binSz

thetas = linspace(binSz/2, 2*pi-binSz/2, length(dirMap))'; % binned angles
rhos   = dirMap ./ nanmax(dirMap(:));   % normalised rates

% get mean vector
rVect = nansum(rhos .* exp(1i*thetas)); % rayleigh vector, complex
meanR = abs(rVect)/nansum(rhos);        % magnitude (real part normalised by weights) 
meanDir = angle(rVect);                 % angle = imaginary part

end

