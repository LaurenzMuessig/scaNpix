function dataSmoothed = smoothWithNaNs(data,kernel)
% smoothWithNaNs - smooth data (e.g. positions) while taking care of NaNs
% in it. Using e.g. 'im_filter' will cause NaNs to spread in the data
% package: scanpix.helpers
%
% Syntax: dataSmoothed = scanpix.helpers.smoothWithNaNs(data,kernel)
%         
%
% Inputs:
%    data   - data that you like to smooth. Must be a column vector and nDims = max 2
%    kernel - smoothing kernel. Must be a column vector              
%
% Outputs:
%    dataSmoothed - your smoothed data
%
% LM 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    data {mustBeNumeric}
    kernel (:,1) {mustBeNumeric}
end

%%
nDims             = size(data,2);

% find NaNs
nanInd            = any(isnan(data),2);

% replace NaNs and do conv
data(nanInd,:)    = 0;
nanFact           = ones(size(data));
nanFact(nanInd,:) = 0;

% convolve data with kernel and get a normalisiation factor for each sample based on n of NaNs in each window
if nDims == 1
    tmpSmooth     = conv(data,kernel,'same');
    nanFact       = conv(nanFact,kernel,'same');
elseif nDims == 2
    tmpSmooth     = conv2(data,kernel,'same');
    nanFact       = conv2(nanFact,kernel,'same');
else
    error(['scaNpix:helpers::smoothWithNaNs: You can have max 2 dimensions in your data - yours has' num2str(nDims)]);
end

% smoothed data
dataSmoothed           = tmpSmooth ./ nanFact;
dataSmoothed(nanInd,:) = NaN; % reassign NaNs back

end