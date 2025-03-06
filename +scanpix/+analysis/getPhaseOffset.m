function [phaseOffset, dirOffset] = getPhaseOffset(crossCorr,corrThr,debugFlag)
% Calculate the offset in spatial phase between 2 grid cells. According to
% Tocker et al., Hippocampus (2015)
% In practice works by Voronoi segementing the spatial autocorr. of the
% spatial crosscorr. and then finding the peak shift in the region of the
% central Voronoi cell in the spatial crosscorr.
%
% package: scanpix.analysis
%
%
% Syntax:
%       phaseOffset = getPhaseOffset(crossCorr)
%       phaseOffset = getPhaseOffset(crossCorr,corrThr)
%       phaseOffset = getPhaseOffset(crossCorr,corrThr,debugOn)
%
% Inputs:
%    crossCorr   - spatial cross corr of 2 grid cells
%    corrThr     - Peaks are r > this value (default=0.1)
%    debugOn     - true/false - plot figure with results (default=false)
%
% Outputs:
%   phaseOffset  - phase offset ([0 pi])
%   dirOffset    - direction of the offset with respect to centre ([0 2pi])
%
%
% LM 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    crossCorr {mustBeNumeric}
    corrThr (1,1) {mustBeNumeric} = 0.1; 
    debugFlag (1,1) {mustBeNumericOrLogical} = false;  
end

%% Voronoi segmentation of autocorr from cross corr
% get auto corr of cross corr
ccAuto = scanpix.analysis.spatialCrosscorr(crossCorr, crossCorr);

% only use central portion of auto corr
centralPoint = ceil([size(crossCorr,2)/2,size(crossCorr,1)/2]); %m,n pair
dx           = -centralPoint(1):centralPoint(2);
RC0          = ceil(size(ccAuto)/2);
ccAuto       = ccAuto(RC0(1)+dx, RC0(2)+dx);
% use watershed to segment AC into peaks
fieldsLabel  = watershed(-ccAuto); %p
stats        = regionprops(fieldsLabel,ccAuto, 'Centroid');
xyCoordMaxBin = reshape([stats.Centroid], 2,[])'; %
% Convert to a new reference frame which as the origin at the centre of the autocorr
xyCoordMaxBinCentral = xyCoordMaxBin-repmat(fliplr(centralPoint), [length(xyCoordMaxBin), 1]);
%find central peak
distFromCentre    = sum(xyCoordMaxBinCentral.^2,2).^0.5;
[~, orderOfClose] = sort(distFromCentre);
closestInd        = orderOfClose(1);

% voronoi segmentation
try
    [V,C] = voronoin(xyCoordMaxBin);
catch
    warning('Can''t compute Voronoi segmentation. Peaks in cross corr are probably not separated well by ''corrThreshold''.');
    [phaseOffset, dirOffset] = deal(NaN);
    return
end

% find peak bin in central Voronoi cell of cross corr
[colsIm, rowsIm]   = meshgrid(1:size(crossCorr,1), 1:size(crossCorr,2));
in                 = inpolygon(rowsIm,colsIm,V(C{closestInd},1),V(C{closestInd},2));
%
crossCorrTemp      = crossCorr;
crossCorrTemp(~in) = 0;
%
[maxVal, maxInd]   = max(crossCorrTemp(:),[],'omitnan');
if maxVal < corrThr
    phaseOffset = pi;
    dirOffset   = NaN;
    warning('No peak found in central Voronoi cell - phaseOffset set to %c.',960);
    return
end
%
[closestPeakCoords(2), closestPeakCoords(1)] = ind2sub(size(crossCorrTemp),maxInd);
% this is the line between centre and AC edge, cutting through 'closestPeakCoords' 
dirOffset = mod(atan2(closestPeakCoords(2)-centralPoint(2),closestPeakCoords(1)-centralPoint(1)),2*pi);
xTmp      = centralPoint(1) * cos(dirOffset) + centralPoint(1);
yTmp      = centralPoint(2) * sin(dirOffset) + centralPoint(2);

% need to loop over all voronoi cell sides
distVoronoi = NaN;
for i = 1:length(C{closestInd})
    
    if i == length(C{closestInd})
        ind2 = 1;
    else
        ind2 = i+1;
    end
    % intersection of voronoi cell and line through peak bin - can break when we found intersection
    [xi,yi] = linexline([centralPoint(1) xTmp], [centralPoint(2) yTmp], [V(C{closestInd}(i),2) V(C{closestInd}(ind2),2)], [V(C{closestInd}(i),1),V(C{closestInd}(ind2),1)]);
    
    if ~isnan(xi); distVoronoi = ( (centralPoint(1)-xi)^2 + (centralPoint(2)-yi)^2 )^0.5; break; end      
end

% phase offset
centreDist  = ( (closestPeakCoords(1)-centralPoint(1))^2 + (closestPeakCoords(2)-centralPoint(2))^2 )^0.5;
phaseOffset = pi*(centreDist/distVoronoi);

%% plot results for checking/debugging
if debugFlag
    figure;
    subplot(1,2,1);
    imagesc(crossCorr); colormap(jet); 
    axis square
    axis off
    hold on
    voronoi(xyCoordMaxBin(:,2),xyCoordMaxBin(:,1));
    plot([centralPoint(1) closestPeakCoords(1)],[centralPoint(2) closestPeakCoords(2)],'k-');
    hold off
    
    subplot(1,2,2);
    imagesc(ccAuto); colormap(jet); 
    axis square
    axis off
    hold on
    voronoi(xyCoordMaxBin(:,2),xyCoordMaxBin(:,1));
    plot([ceil(size(ccAuto,2)/2) xi],[ceil(size(ccAuto,1)/2) yi],'k-');
    hold off
 
end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xi,yi] = linexline(L1x, L1y, L2x, L2y)
%%***********************************************************************%
%*                    Line to line interection point                    *%
%*            Finds the intersection of the two line segments.           *%
%*                                                                      *%
%*                                                                      *%
%* Author: Preetham Manjunatha                                          *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 02/08/2022                                                     *%
%************************************************************************%
%
%************************************************************************%
%
% Usage: [xi,yi] = linexline(L1x, L1y, L2x, L2y)
%
% Inputs:
%
%           L1x                     - Line 1 x1 and x2 coordinates [x1, x2]
%           L1y                     - Line 1 y1 and y2 coordinates [y1, y2]
%           L2x                     - Line 2 x3 and x4 coordinates [x3, x4]
%           L2y                     - Line 2 y3 and y4 coordinates [y3, y4]
%           showIntersectionPlot    - show intersection plot (0 or 1)
% 
% Outputs: 
%
%           xi          - Interection point, x coordinate (NaN if no
%                         interesction)
%           yi          - Interection point, y coordinate (NaN if no
%                         interesction)
%--------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------
% nargin check
if nargin < 4
    error('Not enough input arguments.');
elseif nargin > 5
    error('Too many input arguments.');
end

%------------------------------------------------------------------------------------------------------------------------
% Data
x1 = L1x(1);
y1 = L1y(1);
x2 = L1x(2);
y2 = L1y(2);
x3 = L2x(1);
y3 = L2y(1);
x4 = L2x(2);
y4 = L2y(2);
%------------------------------------------------------------------------------------------------------------------------
% Line segments intersect parameters
u = ((x1-x3)*(y1-y2) - (y1-y3)*(x1-x2)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
%------------------------------------------------------------------------------------------------------------------------
% Check if intersection exists, if so then store the value
if (u >= 0 && u <= 1.0) && (t >= 0 && t <= 1.0)
    xi = ((x3 + u * (x4-x3)) + (x1 + t * (x2-x1))) / 2; 
    yi = ((y3 + u * (y4-y3)) + (y1 + t * (y2-y1))) / 2;
else
    xi = NaN;
    yi = NaN;
end

end


% numerator = abs((V(C{closestInd}(ind2),1) - V(C{closestInd}(i),1)) * (V(C{closestInd}(i),2) - centralPoint(2)) - (V(C{closestInd}(i),1) - centralPoint(1)) * (V(C{closestInd}(ind2),2) - V(C{closestInd}(i),2)));
%         % Find the denominator for point-to-line distance formula.
%         denominator = sqrt((V(C{closestInd}(ind2),1) - V(C{closestInd}(i),1)) ^ 2 + (V(C{closestInd}(ind2),2) - V(C{closestInd}(i),2)) ^ 2);
%         % Compute the distance.
%         distVoronoi = numerator ./ denominator; 

