function IUratio = computeIURatio(ephysObj,trialInd,ellipseFits,cellFilter)
% computeIURatio - compute Intersection/Union ratio inner ring of spatial
% autocorrelation
% package: scanpix.analysis
%

%
% Syntax:
%       IUratio = scanpix.analysis.computeIURatio(obj,trialN)
%       IUratio = scanpix.analysis.computeIURatio(obj,trialN,ellipseFits)
%       IUratio = scanpix.analysis.computeIURatio(obj,trialN,cellFilter)
%
% Inputs:
%    ephysObj    - scanpix.ephys object
%    trialInd    - numeric index for trial which you want to use
%    ellipseFits - array of ellipse fits (one for each cell); should be a nCell by 3 array [orientation, length major axis, length minor axis] (optional)  
%    cellFilter  - logical index for cells you want to use/skip (optional)
%
%
% Outputs: 
%
%    IUratio     - Intersection/union ratio for cell pairs (nCells by nCells)
%
% LM 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% some pre checking

if isempty(ephysObj.maps.sACs{trialInd})
    scanpix.maps.addMaps(ephysObj,'sac',trialInd,ephysObj.mapParams);
end

if nargin < 3 || isempty(ellipseFits)
    [~,~,propsEllipse] = cellfun(@(x) scanpix.analysis.gridprops(x,'getellgridness',true),dataObj.maps.sACs{p.Results.trialn},'uni',0);
    ellipseFits        = cell2mat(cellfun(@(x) [x.ellOrient x.ellAbScale],propsEllipse,'uni',0));
end

if nargin < 4
    cellFilter = true(size(ephysObj.cell_ID,1),1);
end


%% compute ratio

% some init vals for the ellipses
[cols, rows] = meshgrid(1:length(ephysObj.maps.sACs{trialInd}{1}), 1:length(ephysObj.maps.sACs{trialInd}{1}));
center       = ceil(length(ephysObj.maps.sACs{trialInd}{1})/2);

% Intersection/Union ratio - according to Tocker et al. (2015)
IUratio = nan(length(ephysObj.maps.sACs{trialInd}),length(ephysObj.maps.sACs{trialInd}));
for i = 1:length(ephysObj.maps.sACs{trialInd})
    
    if ~cellFilter(i) || any(isnan(ellipseFits(i,2:3))); continue; end
    % Create a logical mask of an ellipse with major/minor axis 'ellipseFits(i,2)', 'ellipseFits(i,3)' and tilt ellipseFits(i,1)
    radiusX = ellipseFits(i,3);
    radiusY = ellipseFits(i,2);
    ellipsePixA = (sin(ellipseFits(i,1)).*(cols - center) + cos(ellipseFits(i,1)).*(rows - center)).^2 ./ radiusX^2 + (cos(ellipseFits(i,1)).*(cols - center) - sin(ellipseFits(i,1)).*(rows - center)) .^2 ./ radiusY^2 <= 1;
    
    for j = i+1:length(ephysObj.maps.sACs{trialInd})
        
        if ~cellFilter(j) || any(isnan(ellipseFits(j,2:3))); continue; end
        
        radiusX = ellipseFits(j,2);
        radiusY = ellipseFits(j,3);
        ellipsePixB = (sin(ellipseFits(j,1)).*(cols - center) + cos(ellipseFits(j,1)).*(rows - center)).^2 ./ radiusX^2 + (cos(ellipseFits(j,1)).*(cols - center) - sin(ellipseFits(j,1)).*(rows - center)) .^2 ./ radiusY^2 <= 1;
        
        IUratio(i,j) = sum(ellipsePixA(:) & ellipsePixB(:)) / sum(ellipsePixA(:) | ellipsePixB(:)); % Intersection / Union
    end
end


end

