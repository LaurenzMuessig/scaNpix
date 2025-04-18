function [bits_per_spike, bits_per_sec] = spatial_info(rMaps, pMap)
% spatial_info - Find spatial information (skaggs info) of place field
% package: scanpix.analysis
%
% Returns Skaggs et al's estimate of spatial information in bits per second:
%
%   I = sum_x p(x) r(x) log(r(x)/r)  
%
% then divide by mean rate over bins to get bits per spike.
%
% Binning could be over any single spatial variable (e.g. location, direction, speed).
% 
% This is a vectorised version of the original fnct in SCAN by TW which increases speed
% and also removes necessity for loop to get spat. info for a bunch of maps
%
%  Usage:  
%    [bits_per_spike, bits_per_sec] = scanpix.analysis.spatial_info(rMaps, pMaps)
%
%  Inputs: 
%    rMaps - single rate map or cell array of rate maps
%    pMaps - single corresponding position map or cell array of corresponding position maps
%
%  Output: 
%    bits_per_spike - 1:nMaps  array of spatial info per spike
%    bits_per_sec   - 1:nMaps  array of spatial info per second
%
% LM 2020
%%
arguments
    rMaps
    pMap {mustBeNumeric}
end

%% format inputs
if ~iscell(rMaps)
    rMaps = {rMaps};
end

% needs to be column vector
% if size(meanRates,2) > size(meanRates,1)
%     meanRates = meanRates';
% end

%% grab a few things
nCells         = length(rMaps); % number of cells
mapSz          = numel(rMaps{1}); % n bins in rate map
duration       = sum(pMap(:), 'omitnan'); % duration of trial

%% get spatial info
rates          = reshape([rMaps{:}],[mapSz nCells]); % rate maps as 3D map x cell arrays
% rates(isnan(rates)) = NaN;
pos            = pMap(:) .* ones(size(rates));
% pos(isnan(pos)) = NaN;
% calculate bits for formula
mean_rate      = sum(rates.*pos,1,'omitnan') ./ duration;
p_x            = pos ./ duration;
% p_r            = bsxfun(@rdivide, rates, meanRates'); 
p_r            = bsxfun(@rdivide, rates, mean_rate); 
%
bits_per_sec   = sum(p_x .* rates .* log2(p_r), 1,'omitnan')';   % sum( p_pos .* rates .* log2(p_rates) )
% bits_per_spike = bits_per_sec ./ meanRates;
bits_per_spike = bits_per_sec ./ mean_rate';

end




