function [meanR,meanDir,thetas,rhos] = rayleighVect(dirMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
binSz  = 360 / length(dirMap); % binSz in degrees;

thetas = linspace(binSz/2, 360-binSz/2, length(dirMap))' .* pi/180; % binned angles
rhos   = dirMap ./ nanmax(dirMap(:)); % normalised rates

% get mean vector
rVect = nansum(rhos .* exp(1i*thetas)); % all vectors, complex
meanR = abs(rVect)/nansum(rhos); %real part normalised by weights 
meanDir = angle(rVect); %complex part=angle

end

