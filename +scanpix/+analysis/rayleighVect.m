function [meanR,meanDir,thetas,rhos] = rayleighVect(dirMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
binSz  = 2*pi / length(dirMap); % binSz in rad;

thetas = linspace(binSz/2, 2*pi-binSz/2, length(dirMap))'; % binned angles
rhos   = dirMap ./ nanmax(dirMap(:));   % normalised rates

% get mean vector
rVect = nansum(rhos .* exp(1i*thetas)); % rayleigh vector, complex
meanR = abs(rVect)/nansum(rhos);        % magnitude (real part normalised by weights) 
meanDir = angle(rVect);                 % angle = imaginary part

end

