function [gridness, Props] = gridprops(autoCorr,fitEllipse,options)
% Using autocorrelogram calculates grid wavelength, orientation, gridness and x,y
% offset of peaks
% package: scanpix.analysis
%
%       [gridness] = scanpix.analysis.gridprops( autoCorr );
%       [gridness, Props] = scanpix.analysis.gridprops( autoCorr, 'paramName', 'paramValue', .. );
%
% Produces key metrics about a grid cell derived from the spatial autocorrelogram: wavelength,
% orientation, gridness, the mean x,y offset of the six central peaks (if fewer than six peaks are
% found then that number is used to calculate the values). Also returns the coordinates of the peaks
% used to calculate these values [x,y] pairs with origin at centre
%
% Optional input parameters: 
%
%   'corrThr'        - Peaks are r > this value
%   'getprops'       - calculate all grid props; if 0 only gridness is calculated
%   'getellgridness' - get grid props by regularising autoCorr from ellipse
%   'minor'          - 
%   'plot'           - make a nice plot showing all grid properties
%   'ax'             - axis handle in case you want plot somewhere specific

%
%  Fields of additional output properties structure ('Props'):
%
%  .gridness
%  .waveLength       NB. Unit for wavelength is bins of autocorr
%  .waveLengthFull   NB. Unit for wavelength is bins of autocorr
%  .orientation
%  .orientationFull
%  .offSet         
%  .fieldSize
%  .closestPeaksCoord 
%
% Note: In the autocorrelogram it is sensible to exclude bins that were constructed with relatively
% small overlap between the ratemap1 and ratemap2 (Hafting excludes bins with an overlap of 20 or
% less). Set these bins to 0 before passing to this function
%
% this function is a re-write of the original gridness calculation written by Tom and Caswell which for several properties didn't work very well in case 
% the grid pattern was slightly irregular. I have tried to more or less replaicate what the Mosers are using.
% 
% LM 2022

%% Params
arguments
    autoCorr {mustBeNumeric}
    fitEllipse (1,1) {mustBeNumericOrLogical} = false;
    options.thresh (1,1) {mustBeNumeric} = -1;
    options.binAC (1,1) {mustBeNumericOrLogical} = false;
    options.nBinSteps (1,1) {mustBeNumeric} = 20;
    options.minPeakSz (1,1) {mustBeNumeric} = 4;
    options.plotEllipse (1,1) {mustBeNumericOrLogical} = false;
    options.ax  {ishghandle(options.ax, 'axes')}
    options.verbose (1,1) {mustBeNumericOrLogical} = false;
    options.legacyMode (1,1) {mustBeNumericOrLogical} = false;
end

%% 
if options.legacyMode
    [gridness, Props] = scanpix.analysis.gridprops_legacy(autoCorr);
    return
end

%%
gridness                = NaN;
%
Props.gridness          = NaN;
Props.wavelength        = NaN;
Props.wavelengthFull    = nan(3,1);
Props.orientation       = NaN;
Props.orientationFull   = nan(3,1);
Props.offset            = NaN;
Props.offsetFull        = nan(3,1);
Props.fieldSize         = NaN;
Props.closestPeaksCoord = nan(6,2);
Props.centralPeakMask   = nan(size(autoCorr));
Props.gridMask          = nan(size(autoCorr));
Props.peakMask          = nan(size(autoCorr));
Props.acReg             = nan(size(autoCorr));
Props.isEllipseFit      = fitEllipse;
Props.ellOrient         = NaN;
Props.ellAbScale        = nan(1,2);
% tmpProps.ax_ind            = nan(3,1);

%%
% Protect against all(autocorr==nan) (when rate map is 0Hz) %
if all(isnan(autoCorr));   return;     end

%%
% --------------------------------------------------------------------------------------------------
% ---- FIND PEAKS ----------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------

% ---- Find central peak in auto corr --- %
% autoCorrTemp = autoCorr;
[xyCoordMaxBin, xyCoordMaxBinCentral, distFromCentre,peakStats, peakMask] = findGridPeaks(autoCorr, options.thresh, options.binAC, options.nBinSteps, options.minPeakSz); 

if isempty(peakStats) || peakStats(1).MajorAxisLength > length(autoCorr)
    if options.verbose; warning('scaNpix::analysis::gridprops: No peaks found. Skipping grid cell properties calculation'); end
    return
end

% Regularise by fitting ellipse to AC 
if fitEllipse
    if options.plotEllipse
        if ~isfield(options,'ax'); options.ax = axes; end
    else
        options.ax = [];
    end
%
    [ ~, ~, orient, abScale ]      = gridEllipse_fit( autoCorr, xyCoordMaxBin(2:end,:), options.plotEllipse,  options.ax, options.verbose );
    if ~isnan(abScale)
        autoCorr                   = regularise_eliptic_grid( autoCorr, abScale, orient*180/pi  );
        % update peak positions
        [xyCoordMaxBin, xyCoordMaxBinCentral, distFromCentre,peakStats, peakMask] = findGridPeaks(autoCorr,options.thresh,  options.binAC, options.nBinSteps, options.minPeakSz);
        if isempty(peakStats) || peakStats(1).MajorAxisLength > length(autoCorr)
            if options.verbose; warning('scaNpix::analysis::gridprops: No peaks found. Skipping grid cell properties calculation'); end
            return
        end
        %
        Props.ellOrient            = orient;
        Props.ellAbScale           = abScale;
        Props.acReg                = autoCorr;
    else
        if options.verbose; warning('scaNpix::analysis::gridprops: Couldn''t fit am ellipse to peaks of autocorr'); end
        return;
    end
end
 
% make central peak mask
[colsIm, rowsIm]                            = meshgrid(1:size(autoCorr,2), 1:size(autoCorr,1));
distMap                                     = sqrt((rowsIm-xyCoordMaxBin(1,1)).^2 + (colsIm-xyCoordMaxBin(1,2)).^2);
centrPeakMask                               = distMap < peakStats(1).MajorAxisLength/2;

% make annulus mask
if ~isempty(distFromCentre)
    maxDist                = max(distFromCentre+mean([peakStats.MajorAxisLength])'/2); 
else
    maxDist                = floor(length(autoCorr)/2); % set to full AC radius in case no peaks are found
end
annMask                    = ~centrPeakMask;
annMask(distMap > maxDist) = false;
%
Props.centralPeakMask      = centrPeakMask;
Props.peakMask             = peakMask | centrPeakMask;
Props.gridMask             = annMask;

%%
% --------------------------------------------------------------------------------------------------
% ---- GRIDNESS ------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
rotAngle = [60, 120, 30, 90, 150];
annCorr  = nan(1,length(rotAngle));
% loop over rotations
for i=1:length(rotAngle)
    autoCorr_rot = imrotate(autoCorr, rotAngle(i), 'bilinear', 'crop');
    nanMask      = ~isnan(autoCorr) & ~isnan(autoCorr_rot); % 
    annCorr(i)   = corr2(autoCorr(annMask & nanMask),autoCorr_rot(annMask & nanMask));
end

if any(isnan(annCorr))
    if options.verbose; warning('scaNpix::analysis::gridprops: Some of the correlations for computing gridness return NaNs - likely means that the central peak is too big to return any meaningful values.'); end
    return;
end

% gridnesss (as per ususal)
gridness       = min(annCorr(1:2)) - max(annCorr(3:end));
Props.gridness = gridness;

%%
% --------------------------------------------------------------------------------------------------
% ---- GRID PROPERTIES -----------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
if length(xyCoordMaxBin) == 7

    % orientation - define 3 axes Moser style
    orientation           = atan2(xyCoordMaxBinCentral(2:end,2),xyCoordMaxBinCentral(2:end,1)); %

    [~, sortInd]          = sort(abs(orientation));
    orientation           = -orientation([sortInd(1:3)]);  % need to flip orientation due to inverted y-axis

    Props.orientation     = circ_mean(orientation);
    Props.orientationFull = orientation;
    % Props.ax_ind          = sortInd(1:3);


    % wavelength
    Props.wavelength      = mean(distFromCentre([sortInd(1:3)]),'omitnan'); % wavelength
    Props.wavelengthFull  = distFromCentre([sortInd(1:3)]);

    % offset
    Props.offsetFull      = pi/4 - abs(mod(orientation(1:min([length(orientation), 3])),pi/2) - pi/4); % from Stensola et al. (2015)
    Props.offset          = min( Props.offsetFull ); %
    %
    tmp                     = ceil(xyCoordMaxBin(2:end,:));
    Props.closestPeaksCoord = tmp(sortInd,:);

    % field size THIS NEEDS WORK AS NOT MATCHING WITH THE GENERAL ALGORTIHM
    Props.fieldSize = pi*(peakStats(1).EquivDiameter/2)^2;   
else
    if options.verbose; warning('scaNpix::analysis::gridprops:Not enough peaks detected for grid property calculation.'); end
end



end

% -------------------------------------------------------------------------------------------------
% --- INLINE FUNCTIONS ----------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
function [xyCoordMaxBin, xyCoordMaxBinCentral, distFromCentre, peakStats, peakMask] = findGridPeaks(autoCorr,thresh,binAC,nBinSteps,minPeakSz)
%%
% autoCorrTemp = autoCorr;

[peakStats, peakMask] = scanpix.analysis.fieldDetect(autoCorr,'thrMode','abs','thr',thresh,'binEdges',[-1 1],'binMap',binAC,'minPeakSz',minPeakSz,'nBinSteps',nBinSteps,'debugOn', false);

%
xyCoordMaxBin        = round(reshape([peakStats.WeightedCentroid], 2,[])'); %Still x,y pair

centralPoint         = ceil([size(autoCorr,2)/2,size(autoCorr,1)/2]); %m,n pair
xyCoordMaxBinCentral = xyCoordMaxBin-repmat(fliplr(centralPoint), [size(xyCoordMaxBin,1), 1]);
distFromCentre       = sum(xyCoordMaxBinCentral.^2,2).^0.5;
% exit gracefully
if ~any(peakMask(:)); return; end

% now get the 6 closest peaks to cental one
[~, orderOfClose]     = sort(distFromCentre);
%
orderOfClose          = orderOfClose(1:min([7,length(orderOfClose)])); % take central peak + 6 closest peaks
distFromCentre        = distFromCentre(orderOfClose);
%
xyCoordMaxBin         = xyCoordMaxBin(orderOfClose,:);
xyCoordMaxBinCentral  = xyCoordMaxBinCentral(orderOfClose,:);
peakStats             = peakStats(orderOfClose);
% remove all but closest 6 peaks from peak mask
allInd                      = 1:numel(autoCorr);
setToZero                   = ~ismember(allInd',vertcat(peakStats.PixelIdxList));
peakMask(allInd(setToZero)) = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From Barry lab:
function [ xyScale, eccent, orient, abScale ] = gridEllipse_fit( sac, closestPeaksCoord, plotEllipse, hAx, verboseFlag )
%GRIDELLIPSE_FIT Fits elipse to grid sac - esimates xy scale
% Grids are often not regular. This code takes a sac and attempts to fit an
% elipse to the central six peaks. This enables an estimate of the scale in
% x and y dimenson as well as measurement of how eliptical the grid is.
% Scale in x dim is defined as point in which the elipse passes through y=0
% and vice versa
%
% IMPORTANT. Code must find 6 peaks in the sac otherwise it cannot fit an
% ellipse using the least squares method (even though one is actually
% defined). In these situations all values are returned as nan.
%
% WARNING. Not fully tested values that this function produces. xyScale
% seems to be broadly correct but might be out by a small factor. Requires
% further testing before publishing results.
% 
%
% ARGS
% sac      spatial autocorr of a grid, best to construct from smoothed ratemap
% showElipseFig [not required - default 'false'] 'true' or 'false'
%
% RETURNS
% xyScale [1x2] - scale in x dimension and scale in y dimension in bins
% eccent [1] - eccentricity of the elipse, 0=circular.
% orient [1] - orientation of major axis anti-cw from x-axis in rads [for
%           ij orientation of image]
% abScale [1x2] - scale of major and minor axis (a is major, b is minor]
%
% EXAMPLE
% [ xyScale, eccent, orient, abScale ] = gridElipse_fit(sac, true, 4);
% [ xyScale, eccent, orient, abScale ] = gridElipse_fit(sac);


% -------------------------------------------------------------------------
% --- PARSE ARGUMENTS AND VARIABLES ---------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% --- MAIN FUNCTION -------------------------------------------------------
% -------------------------------------------------------------------------
%FIRST BLOCK
%Note first blocks of code borrows heavily from autoCorrProps and does
%basic processing of SAC to get six central points not including the
%central peak.
% 
% sac(isnan(sac))=-1; %Sub nans for -1
% autoCorrTemp = sac;        % TW. Do not allow local max that have .. [lines addopted from TW code]
% autoCorrTemp(sac<=0) = -1; %  .. r-value below zero.
% %Don't consider imaginary components which some times appear in shuffled data
% autoCorrTemp=real(autoCorrTemp);
% peaksAutoCorr= imregionalmax(autoCorrTemp); %Find local maxima
% clear autoCorrTemp
% 
% [lableMask, ~]=bwlabel(peaksAutoCorr, 8);
% 
% %In case adjacent points share maxima find the centroid of them - NB returns structure array
% %stats(n).Centroid containing for each peak the x,y position of the centroid but y is
% %counting down from origin at top left
% stats=regionprops(lableMask, 'Centroid');
% %NL. [n x 2] pairs of x,y coord for max points
% xyCoordMaxBin=reshape([stats.Centroid], 2,[])'; %Still x,y pair
% 
% % Convert to a new reference frame which as the origin at the centre of the autocorr
% % and with y negative at top and y positive at bottom, x negative on left and postive on
% % right
% % NB autocorr will always have sides with odd number of bins
% centralPoint=ceil(size(sac)/2); %m,n pair
% xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPoint), [size(xyCoordMaxBin,1), 1]);
% 
% %Calculate distance of peaks from centre point and find seven closest (one will be central peak
% %disregard this)
% distFromCentre=sum(xyCoordMaxBinCentral.^2,2).^0.5;
% [~, orderOfClose]= sort(distFromCentre);
% 
% %Get id of closest peaks - note closest peak 1 will be centre
% if length(orderOfClose)>=7; closestPeaks=orderOfClose(1:7); %Might be fewer than 7 peaks
% else closestPeaks=orderOfClose(1:end);
% end
% 
% %x,y pairs in cartesian coords with y counting down from top and origin at top left.
% closestPeaksCoord=xyCoordMaxBin(closestPeaks,:);
% closestPeaksCoord=closestPeaksCoord([2:end],:); %But not central one

%Check how many peaks are found - must be ==6 to proceed
if length(closestPeaksCoord)<6
    [ xyScale, eccent, orient, abScale ] = deal(nan);
    if verboseFlag; warning('scaNpix::analysis::gridprops:Too few peaks to define elipse - returning nan'); end
    return
end


% SECOND BLOCK
% Fit ellipse to the central points

%Option to draw elipse onto sac - useful for debugging - do this is second
%arg is true
if ~plotEllipse %Don't draw
    elipseData=sf_fit_ellipse(closestPeaksCoord(:,1), closestPeaksCoord(:,2),[],verboseFlag);
elseif plotEllipse %Do draw
    imagesc(hAx,sac); %draw sac
    axis(hAx,'equal','off');
    hold(hAx,'on');
    scatter(hAx,closestPeaksCoord(:,1), closestPeaksCoord(:,2)); %Draw on peaks
    elipseData=sf_fit_ellipse(closestPeaksCoord(:,1), closestPeaksCoord(:,2), hAx, verboseFlag);
    hold(hAx,'off');
    % set(hAx,'ydir','normal');
    title(hAx,'Ellipse fit');
end

%Check if an ellipse was fit - if not elipseData will be empty
if isempty(elipseData) || isempty(elipseData.a) %No ellipse found; LM EDIT!!! 
    if verboseFlag; warning('scaNpix::analysis::gridprops:Failed to fit ellipse - returning nan'); end
    [ xyScale, eccent, orient, abScale ] = deal(nan);
    return
end

%Pull out data about elipse from structure returned by sub functions
a=elipseData.a; %Length of major axis - conventionally 'a'
b=elipseData.b; %Length of minor axis - conventionally 'b'
abScale=[a,b];

%NL is orient of major axis anticlockwise from x-axis in rads but note this
%is for mn coordinates (i.e. origin top left) so in conventional xy
%coordinates there should be a negative sign in front of this value
orient=elipseData.phi; % orient of major axis anticlock from x-axis in rads

%Code sometimes flips a and b
if a<b %Flipped
    b=elipseData.a;
    a=elipseData.b;
    orient=mod(elipseData.phi+(pi/2), 2*pi); %Add 90deg to orient
end

eccent=sqrt((a^2 - b^2)/a^2 ); %Eccentricity where 0 is a circle

%Now get xy scale - use equation for elipse to find value of x when y=0 &
%vice versa
%First working in polar coordinates define where x and y aixs would lie on
%equivalent non-rotated elipse. Again note we're working in mn coords
theta_r_axis = [orient, pi/2+orient]; %pol coord in rads equivalent to x-axis and y-axis
ellipse_x = a*cos( theta_r_axis ); %x value for cart coord
ellipse_y =  b*sin( theta_r_axis );%y value for cart coord
xyScale=sqrt(sum([ellipse_x; ellipse_y].^2,1)); %[xScale, yScale]


end

% From Barry lab (based on fxchange function):
function ellipse_t = sf_fit_ellipse( x,y,axis_handle, verboseFlag )
%
% fit_ellipse - finds the best fit to an ellipse for the given set of points.
%
% Format:   ellipse_t = fit_ellipse( x,y,axis_handle )
%
% Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
%           axis_handle - optional. a handle to an axis, at which the estimated ellipse
%                         will be drawn along with it's axes
%
% Output:   ellipse_t - structure that defines the best fit to an ellipse
%                       a           - sub axis (radius) of the X axis of the non-tilt ellipse
%                       b           - sub axis (radius) of the Y axis of the non-tilt ellipse
%                       phi         - orientation in radians of the ellipse (tilt)
%                       X0          - center at the X axis of the non-tilt ellipse
%                       Y0          - center at the Y axis of the non-tilt ellipse
%                       X0_in       - center at the X axis of the tilted ellipse
%                       Y0_in       - center at the Y axis of the tilted ellipse
%                       long_axis   - size of the long axis of the ellipse
%                       short_axis  - size of the short axis of the ellipse
%                       status      - status of detection of an ellipse
%
% Note:     if an ellipse was not detected (but a parabola or hyperbola), then
%           an empty structure is returned

% =====================================================================================
%                  Ellipse Fit using Least Squares criterion
% =====================================================================================
% We will try to fit the best ellipse to the given measurements. the mathematical
% representation of use will be the CONIC Equation of the Ellipse which is:
%
%    Ellipse = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
%
% The fit-estimation method of use is the Least Squares method (without any weights)
% The estimator is extracted from the following equations:
%
%    g(x,y;A) := a*x^2 + b*x*y + c*y^2 + d*x + e*y = f
%
%    where:
%       A   - is the vector of parameters to be estimated (a,b,c,d,e)
%       x,y - is a single measurement
%
% We will define the cost function to be:
%
%   Cost(A) := (g_c(x_c,y_c;A)-f_c)'*(g_c(x_c,y_c;A)-f_c)
%            = (X*A+f_c)'*(X*A+f_c)
%            = A'*X'*X*A + 2*f_c'*X*A + N*f^2
%
%   where:
%       g_c(x_c,y_c;A) - vector function of ALL the measurements
%                        each element of g_c() is g(x,y;A)
%       X              - a matrix of the form: [x_c.^2, x_c.*y_c, y_c.^2, x_c, y_c ]
%       f_c            - is actually defined as ones(length(f),1)*f
%
% Derivation of the Cost function with respect to the vector of parameters "A" yields:
%
%   A'*X'*X = -f_c'*X = -f*ones(1,length(f_c))*X = -f*sum(X)
%
% Which yields the estimator:
%
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       |  A_least_squares = -f*sum(X)/(X'*X) ->(normalize by -f) = sum(X)/(X'*X)  |
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% (We will normalize the variables by (-f) since "f" is unknown and can be accounted for later on)
%
% NOW, all that is left to do is to extract the parameters from the Conic Equation.
% We will deal the vector A into the variables: (A,B,C,D,E) and assume F = -1;
%
%    Recall the conic representation of an ellipse:
%
%       A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
%
% We will check if the ellipse has a tilt (=orientation). The orientation is present
% if the coefficient of the term "x*y" is not zero. If so, we first need to remove the
% tilt of the ellipse.
%
% If the parameter "B" is not equal to zero, then we have an orientation (tilt) to the ellipse.
% we will remove the tilt of the ellipse so as to remain with a conic representation of an
% ellipse without a tilt, for which the math is more simple:
%
% Non tilt conic rep.:  A`*x^2 + C`*y^2 + D`*x + E`*y + F` = 0
%
% We will remove the orientation using the following substitution:
%
%   Replace x with cx+sy and y with -sx+cy such that the conic representation is:
%
%   A(cx+sy)^2 + B(cx+sy)(-sx+cy) + C(-sx+cy)^2 + D(cx+sy) + E(-sx+cy) + F = 0
%
%   where:      c = cos(phi)    ,   s = sin(phi)
%
%   and simplify...
%
%       x^2(A*c^2 - Bcs + Cs^2) + xy(2A*cs +(c^2-s^2)B -2Ccs) + ...
%           y^2(As^2 + Bcs + Cc^2) + x(Dc-Es) + y(Ds+Ec) + F = 0
%
%   The orientation is easily found by the condition of (B_new=0) which results in:
%
%   2A*cs +(c^2-s^2)B -2Ccs = 0  ==> phi = 1/2 * atan( b/(c-a) )
%
%   Now the constants   c=cos(phi)  and  s=sin(phi)  can be found, and from them
%   all the other constants A`,C`,D`,E` can be found.
%
%   A` = A*c^2 - B*c*s + C*s^2                  D` = D*c-E*s
%   B` = 2*A*c*s +(c^2-s^2)*B -2*C*c*s = 0      E` = D*s+E*c
%   C` = A*s^2 + B*c*s + C*c^2
%
% Next, we want the representation of the non-tilted ellipse to be as:
%
%       Ellipse = ( (X-X0)/a )^2 + ( (Y-Y0)/b )^2 = 1
%
%       where:  (X0,Y0) is the center of the ellipse
%               a,b     are the ellipse "radiuses" (or sub-axis)
%
% Using a square completion method we will define:
%
%       F`` = -F` + (D`^2)/(4*A`) + (E`^2)/(4*C`)
%
%       Such that:    a`*(X-X0)^2 = A`(X^2 + X*D`/A` + (D`/(2*A`))^2 )
%                     c`*(Y-Y0)^2 = C`(Y^2 + Y*E`/C` + (E`/(2*C`))^2 )
%
%       which yields the transformations:
%
%           X0  =   -D`/(2*A`)
%           Y0  =   -E`/(2*C`)
%           a   =   sqrt( abs( F``/A` ) )
%           b   =   sqrt( abs( F``/C` ) )
%
% And finally we can define the remaining parameters:
%
%   long_axis   = 2 * max( a,b )
%   short_axis  = 2 * min( a,b )
%   Orientation = phi
%
%

% initialize
orientation_tolerance = 1e-3;

% empty warning stack
warning( '' );

% prepare vectors, must be column vectors
x = x(:);
y = y(:);

% remove bias of the ellipse - to make matrix inversion more accurate. (will be added later on).
mean_x = mean(x);
mean_y = mean(y);
x = x-mean_x;
y = y-mean_y;

% the estimation for the conic equation of the ellipse
X = [x.^2, x.*y, y.^2, x, y ];
a = sum(X)/(X'*X);

% check for warnings
if ~isempty( lastwarn )
    if verboseFlag; warning('scaNpix::analysis::gridprops: Stopped because of a warning regarding matrix inversion' ); end
    ellipse_t = [];
    return
end

% extract parameters from the conic equation
[a,b,c,d,e] = deal( a(1),a(2),a(3),a(4),a(5) );

% remove the orientation from the ellipse
if ( min(abs(b/a),abs(b/c)) > orientation_tolerance )
    
    orientation_rad = 1/2 * atan( b/(c-a) );
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    [a,b,c,d,e] = deal(...
        a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
        0,...
        a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
        d*cos_phi - e*sin_phi,...
        d*sin_phi + e*cos_phi );
    [mean_x,mean_y] = deal( ...
        cos_phi*mean_x - sin_phi*mean_y,...
        sin_phi*mean_x + cos_phi*mean_y );
else
    orientation_rad = 0;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
end

% check if conic equation represents an ellipse
test = a*c;
switch (1)
    case (test>0),  status = '';
    case (test==0), status = 'Parabola found';  if verboseFlag; warning( 'fit_ellipse: Did not locate an ellipse' ); end
    case (test<0),  status = 'Hyperbola found'; if verboseFlag; warning( 'fit_ellipse: Did not locate an ellipse' ); end
end

% if we found an ellipse return it's data
if (test>0)
    
    % make sure coefficients are positive as required
    if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end
    
    % final ellipse parameters
    X0          = mean_x - d/2/a;
    Y0          = mean_y - e/2/c;
    F           = 1 + (d^2)/(4*a) + (e^2)/(4*c);
    [a,b]       = deal( sqrt( F/a ),sqrt( F/c ) );
    long_axis   = 2*max(a,b);
    short_axis  = 2*min(a,b);
    
    % rotate the axes backwards to find the center point of the original TILTED ellipse
    R           = [ cos_phi sin_phi; -sin_phi cos_phi ];
    P_in        = R * [X0;Y0];
    X0_in       = P_in(1);
    Y0_in       = P_in(2);
    
    % pack ellipse into a structure
    ellipse_t = struct( ...
        'a',a,...
        'b',b,...
        'phi',orientation_rad,...
        'X0',X0,...
        'Y0',Y0,...
        'X0_in',X0_in,...
        'Y0_in',Y0_in,...
        'long_axis',long_axis,...
        'short_axis',short_axis,...
        'status','' );
else
    % report an empty structure
    ellipse_t = struct( ...
        'a',[],...
        'b',[],...
        'phi',[],...
        'X0',[],...
        'Y0',[],...
        'X0_in',[],...
        'Y0_in',[],...
        'long_axis',[],...
        'short_axis',[],...
        'status',status );
end

% check if we need to plot an ellipse with it's axes.
if (nargin>2) && ~isempty( axis_handle ) && (test>0)
    
    % rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    
    % the axes
    ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
    horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    
    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = X0 + a*cos( theta_r );
    ellipse_y_r     = Y0 + b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    % draw
%     hold_state = get( axis_handle,'NextPlot' );
%     set( axis_handle,'NextPlot','add' );
    plot(axis_handle, new_ver_line(1,:),new_ver_line(2,:),'r' );
    plot(axis_handle, new_horz_line(1,:),new_horz_line(2,:),'r' );
    plot(axis_handle, rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
%     set( axis_handle,'NextPlot',hold_state );
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From Barry lab:
function [ regSac ] = regularise_eliptic_grid( sac, abScale, orient  )
%REGULARISE_ELIPTIC_GRID Reshapes eliptical grid to be regular
% Similar to analyses run by other groups - having found elipticalness of
% grid using grid_xy_elipse_scale.m can then use this function to make grid
% points lie on circle and then redo gridness.
%
% Method is to first rotate the sac so that the major axis is aligned to
% the y-axis. Then use the major minor axis to adjust scale. Finally rotate
% back
%
% ARGS
% sac - spatial autocorr
% abScale - scale of major and minor axis
% orient - orientation in deg of major axis antiC from x-axis

%--- House keeping
abScale=sort(abScale, 'descend'); %Sometimes ab scale is mixed up - major shoudl be first


%--- Main function
%First check if need to do anything
if abScale(1)==abScale(2) %Grid is already regular
     regSac=sac;
    return
end

% LM edit
tmpSAC             = sac;
tmpSAC(isnan(sac)) = 0; % leaving the nan's in sac makes them spread during the interpolation

%Grids aren't regular to start to regularise
%1)First rotate so major axis aligns to x-axis
% regSac=imrotate(sac, -orient, 'bilinear');
regSac=imrotate(tmpSAC, -orient, 'bilinear');


%2)Decide by how much to resize - major axis is x so work on y axis to
%bring to same scale
sacSize=size(regSac);
xStart=1;
xEnd=sacSize(2);
tmp=sacSize(1)/2 - (sacSize(1)/2)*(abScale(2)/abScale(1));
yStart=1+tmp;
yEnd=sacSize(1)-tmp;
% LM EDIT: need to add +1 to size in case it's even
sacSize(rem(sacSize,2) == 0) = sacSize(rem(sacSize,2) == 0) + 1;    
%3) Do the resample the sac to regularise the grid
% regSac=interp2(regSac, linspace(xStart,xEnd, sacSize(2)), linspace(yStart, yEnd, sacSize(1))');
% LM edit
regSac=interp2(regSac, linspace(xStart,xEnd, sacSize(2)), linspace(yStart, yEnd, sacSize(1))'); 


%4) Rotate back to original orientation
% regSac=imrotate(regSac, orient, 'bilinear','crop');
% LM edit - crop to keep size == odd
if size(tmpSAC,1) == size(tmpSAC,2)
    regSac=imrotate(regSac, orient, 'bilinear','crop');
else
    % in case ac isn't square we can't use crop
    regSac=imrotate(regSac, orient, 'bilinear');
end

%5) Finally keep only the central portion of the sac so that it matches
%the original size
sizeDif=round((size(regSac)-size(sac))/2);
regSac=regSac(1+sizeDif(1):end-sizeDif(1), 1+sizeDif(2):end-sizeDif(2));
% LM edit
regSac(regSac==0) = NaN; % reassign original nan's 

end


    
    % don't allow peaks which are too close - only keep closer one to centre.
%     [r,c] = find(triu( abs(circ_dist2(orientation)) < minOr, 1 ));
%     delInd = [];
%     for i = 1:length(r)
%         if distFromCentre(r(i)) < distFromCentre(c(i))
%             delInd = [delInd;c(i)];
%         else
%             delInd = [delInd;r(i)];
%         end
%     end
%     delInd = unique(delInd);
%     % remove peaks
%     orientation(delInd)     = [];
%     xyCoordMaxBin(delInd,:) = [];
%     distFromCentre(delInd)  = [];
%     [~, orderOfClose]       = sort(distFromCentre);
%     orderOfClose            = orderOfClose(1:min([6,length(orientation)])); % in case there >6 peaks keep 6 closest ones
%     % order by distance
%     orientation   = mod(orientation(orderOfClose),2*pi);
%     wavelength    = distFromCentre(orderOfClose); % wavelength
%     xyCoordMaxBin = xyCoordMaxBin(orderOfClose,:);


    % this is from the Moser lab code base
%     orientation = mod(orientationRaw,60*pi/180);
%     orientation(orientation-60*pi/180<abs(orientation)) = orientation-60*pi/180;
%     if circ_mean(abs(abs(orientation)-30*pi/180)) < circ_mean(abs(orientation))
%         if circ_median(orientation) < 0
%             orientation(orientation>0) = -1*orientation(orientation>0);
%         else
%             orientation(orientation<0) = -1*orientation(orientation<0);
%         end
%     end
%     Props.orientation = circ_mean(orientation);
     
%     if min([abs(o-60*pi/180);abs(o)]) == 0
        

%     [~,ax1_ind]  = min(abs(orientation)); 
%     tmp(ax1_ind) = NaN;
%     
%     [~, sortInd] = sort(circ_dist(abs(tmp),abs(orientation(ax1_ind))));
% %     [~, ax2_ind] = min(circ_dist(abs(tmp),abs(orientation(ax1_ind))));
% %     tmp(ax2_ind) = NaN;
% %     [~, ax3_ind] = min(circ_dist(abs(tmp),abs(orientation(ax1_ind))));
%     orientation = orientation([ax1_ind sortInd(1) sortInd(2)]);


% if 0
%     % % this works pretty well
%     % autoCorrTemp(~annMask) = NaN;
%     % zAutoCorr              = (autoCorrTemp - mean(autoCorrTemp(:),'omitnan')) ./ std(autoCorrTemp(:),'omitnan');
%     % peaksAutoCorr          = zAutoCorr > thr;
% 
% else
%     autoCorrTemp(isnan(autoCorrTemp)) = 0;
%     fieldsLabel  = watershed(-autoCorrTemp);
% 
%     thresholds = nan(size(autoCorrTemp));
%     for i = 1:max(fieldsLabel(:))
%         thresholds(fieldsLabel == i) = max(autoCorrTemp(fieldsLabel == i),[],'omitnan')/2;
%     end
%     peaksAutoCorr          = autoCorrTemp > thresholds;
%     peaksAutoCorr = imclose(peaksAutoCorr,strel('square',3));
% end

