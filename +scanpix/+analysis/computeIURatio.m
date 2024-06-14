function IUratio = computeIURatio(spatialACs,ellipseFits)
% computeIURatio - compute Intersection/Union ratio of inner ring of
% spatial autocorrelation across grid cells
% package: scanpix.analysis
%
%
% Syntax:
%       IUratio = scanpix.analysis.computeIURatio(spatialACs)
%       IUratio = scanpix.analysis.computeIURatio(spatialACs,ellipseFits)
%
% Inputs:
%    spatialACs  - cell array of spatial autocorrelations of grid cells
%    ellipseFits - array of ellipse fits (one for each cell); nCell by 3 array [orientation, length major axis, length minor axis] (optional)  
%
% Outputs: 
%
%    IUratio     - Intersection/union ratio for cell pairs (nCells by nCells)
%
% LM 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% some pre-checking
if nargin == 1 
    %%%% NNEDS FIXING
    % [~,~,propsEllipse] = cellfun(@(x) scanpix.analysis.gridprops(x,'getellgridness',true),spatialACs,'uni',0);
    % ellipseFits        = cell2mat(cellfun(@(x) [x.ellOrient x.ellAbScale],propsEllipse,'uni',0));
end

%% compute ratio
% some init vals for the ellipses
[cols, rows] = meshgrid(1:length(spatialACs{1}), 1:length(spatialACs{1}));
center       = ceil(length(spatialACs{1})/2);

% Intersection/Union ratio - according to Tocker et al. (2015)
IUratio = nan(length(spatialACs),length(spatialACs));

for i = 1:length(spatialACs)
    
    if any(isnan(ellipseFits(i,2:3))); continue; end
    % Create a logical mask of an ellipse with major/minor axis 'ellipseFits(i,2)', 'ellipseFits(i,3)' and tilt ellipseFits(i,1)
    radiusX     = ellipseFits(i,3);
    radiusY     = ellipseFits(i,2);
    ellipsePixA = (sin(ellipseFits(i,1)).*(cols - center) + cos(ellipseFits(i,1)).*(rows - center)).^2 ./ radiusX^2 + (cos(ellipseFits(i,1)).*(cols - center) - sin(ellipseFits(i,1)).*(rows - center)) .^2 ./ radiusY^2 <= 1;
    %
    for j = i+1:length(spatialACs)
        
        if any(isnan(ellipseFits(j,2:3))); continue; end
        
        radiusX      = ellipseFits(j,2);
        radiusY      = ellipseFits(j,3);
        ellipsePixB  = (sin(ellipseFits(j,1)).*(cols - center) + cos(ellipseFits(j,1)).*(rows - center)).^2 ./ radiusX^2 + (cos(ellipseFits(j,1)).*(cols - center) - sin(ellipseFits(j,1)).*(rows - center)) .^2 ./ radiusY^2 <= 1;
        %
        IUratio(i,j) = sum(ellipsePixA(:) & ellipsePixB(:),'omitnan') / sum(ellipsePixA(:) | ellipsePixB(:),'omitnan'); % Intersection / Union
    end
end


end

