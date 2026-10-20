function prms = defaultParamsRateMaps
% defaultParamsRateMaps - Generate the default params structure for making rate maps
% package: scanpix.maps
%
% Syntax:  prms = scanpix.maps.defaultParamsRateMaps
%
% Inputs:
%
% Outputs: prms - prms struct with default parameters
%
%
% LM 2020

%% params for maps

%% general
speedFilterLimits                = [2.5 400];
showWaitBar                      = false;

%% rate maps
prms.rate.speedFilterFlagRMaps   = 1;    % y/n
prms.rate.speedFilterLimitLow    = speedFilterLimits(1);
prms.rate.speedFilterLimitHigh   = speedFilterLimits(2);
prms.rate.binSizeSpat            = 2.5;  % spatial bin size in cm^2
prms.rate.smooth                 = 'adaptive'; % 'boxcar; 'adaptive'
prms.rate.kernel                 = 5;    % size smoothing kernel (in bins)
prms.rate.alpha                  = 200;  % alpha parameter - don't change
prms.rate.trimNaNs               = false;
prms.rate.showWaitBar            = showWaitBar;

%% dir maps
prms.dir.speedFilterFlagDMaps    = 1;    % y/n
prms.dir.speedFilterLimitLow     = speedFilterLimits(1);
prms.dir.speedFilterLimitHigh    = speedFilterLimits(2);
prms.dir.binSizeDir              = 6;    % in degrees
prms.dir.dirSmoothKern           = 5;    % smooth across that many bins
prms.dir.showWaitBar             = showWaitBar;

%% lin maps
% prms.lin.trackType               = [];
% Parameters for linearisation of position data 
% prms.linpos.trackLength           = [];   % in pixels; for square track it should be length for 1 arm 
prms.linpos.minDwellForEdge       = 1;    % in s
prms.linpos.durThrCohRun          = 2;   % In seconds
prms.linpos.filtSigmaForRunDir    = 3;    % Units=sec. Sigma of the Gaussian filter to pre-filter the data before finding CW and CCW runs. Kernel is 2*sigma in length.
prms.linpos.durThrJump            = 2;    % In seconds
prms.linpos.gradThrForJumpSmooth  = 2.5;
prms.linpos.runDimension          = [];           % This is the dimension to run (X=1, Y=2); will be estimated if empty (lin track only parameter) 
prms.linpos.dirTolerance          = 70 * pi/180;  % Tolerance for heading direction, relative to arm direction, for calculating run direction (degrees, but will switch to radians)
% linear map params
prms.lin.binSizeLinMaps          = 2.5; % bin size (cm^2)
prms.lin.smoothFlagLinMaps       = 1;  % y/n
prms.lin.smoothKernelSD          = 2;  % SD, in bins
prms.lin.speedFilterFlagLMaps    = 1;  % y/n
prms.lin.speedFilterLimitLow     = speedFilterLimits(1);
prms.lin.speedFilterLimitHigh    = speedFilterLimits(2);
prms.lin.posIsCircular           = 0;
prms.lin.remTrackEnds            = 0;  % Remove this many bins from each end of the linear track. Set to 0 to remove none.
prms.lin.showWaitBar             = showWaitBar;

%% object vector maps
prms.objVect.speedFilterFlagOVMaps = 1;  % y/n
prms.objVect.speedFilterLimitLow   = speedFilterLimits(1);
prms.objVect.speedFilterLimitHigh  = speedFilterLimits(2);
prms.objVect.minDist             = 0;     % in cm; distance bins start here
prms.objVect.maxDist             = [];    % in cm; empty = max distance in data
prms.objVect.binSz_dist          = 2.5;   % in cm;  2cm in Høydal et al (2019)
prms.objVect.binSz_dir           = 5;     % in degrees;  5deg in Høydal et al (2019)
prms.objVect.smKernelSz_OV       = [];    % in bins; empty = 2*ceil(2*smSigma_OV)+1 (as imgaussfilt)
prms.objVect.smSigma_OV          = 2;
prms.objVect.showWaitBar         = showWaitBar;

%% speed maps
prms.speed.minBinProp            = 0.005; % valid speed bins need to contain > prctg of samples of population (0.5%)
prms.speed.binSizeSpeed          = 2;     % cm/s (2 cm/s)
prms.speed.maxSpeed              = 40;    % cm/s (40 cm/s)
% prms.speed.smKernelLength        = 3;     % Note: this is the SD of the kernel!!
% prms.speed.normaliseFR           = true;  % logical flag
prms.speed.confInt               = 95;    % confidence interval
prms.speed.showWaitBar           = showWaitBar;

%% spatial autocorellograms
% using the Mosers way to compute the 2D autocorrelation produces much better results than using the functions that Caswell wrote - both
% methods produce the same AC for some rate maps but not generally. I don't understand why. 
prms.sac.method                  = 'moser'; 
prms.sac.removeMinOverlap        = true;
prms.sac.smooth                  = false;
prms.sac.hSize                   = 5;
prms.sac.sigma                   = 1.5;

%% grid props
% after many rounds of testing and looking at results, it seems the best results for finding the grid peaks in a noisy AC is to use watershedding (with a few tricks to avoid oversegmentation) 
% on a binned version of the AC that is threshholded at r=0 
prms.gridProps.binAC             = true;     % 
prms.gridProps.nBinSteps         = 21;     % bin size = 0.1 
prms.gridProps.thresh            = 0;  % 
prms.gridProps.minPeakSz         = 8;    % 
prms.gridProps.plotEllipse       = false;         
prms.gridProps.verbose           = false;

%% grid props v2 (scanpix.analysis.gridprops_v2, plot: scanpix.plot.plotGridProps_v2)
% Peaks are found via the troughs around the central peak: for a true neighbour peak P the minimum along
% the ray centre->P sits at |P|/2 (any lattice, also elliptical). Candidate peaks failing that test are
% dropped, the remaining ones are combined into the lattice (a1, a2, a1+-a2) with the best peak-vs-hole
% contrast. Separate struct because scanpix.analysis.gridprops errors on unknown options.
% Usage: gpArgs = scanpix.helpers.getGridPropsV2Params(prms); [g, P] = scanpix.analysis.gridprops_v2(sac, gpArgs{:});
% peak finding
prms.gridPropsV2.cRange          = [0.75 1.2]; % allowed |peak| / (2*trough) along a ray; <1 allowed as peaks cut off by the SAC edge read short
prms.gridPropsV2.band            = [1.2 3.2];  % candidate peaks searched within this range x central trough radius
prms.gridPropsV2.minSep          = 35;         % min angle (deg) between two lattice axes
prms.gridPropsV2.wedge           = 3;          % angular average (deg) for the ray profiles
prms.gridPropsV2.edgeTol         = 1;          % peaks within this many bins of the NaN edge of the SAC are flagged 'censored'
% rejection rules (avoid fitting a lattice to noise -> then no peaks/scale/orientation, gridness still computed)
prms.gridPropsV2.minProm         = 0.1;        % peak must exceed the trough on its ray by this (r)
prms.gridPropsV2.minContrast     = 0.25;       % min SAC contrast (predicted peaks - predicted holes); 96% of grid cells (G>0.5) pass, ~half of non-grid cells don't
prms.gridPropsV2.require3        = false;      % true: reject lattice unless all 3 axes have a detected peak (stricter; drops ~9% of grid cells, mostly large scale)
% gridness
prms.gridPropsV2.minAnnWidth     = 1.5;        % first outer annulus radius >= inner radius * (1 + minAnnWidth) (inner radius = trough ring)
prms.gridPropsV2.maxEllRatio     = 1.5;        % ellipse correction always used up to this major/minor ratio ...
prms.gridPropsV2.maxEllRatioSupported = 2.5;   % ... and up to this ratio if all 3 axes have a detected peak not cut off by the SAC edge (else gridnessEll = gridness)
prms.gridPropsV2.ellMode         = 'compress'; % ellipse -> circle by compressing the major axis ('stretch' minor / 'area' both / 'v1' = gridprops style performed the same)
prms.gridPropsV2.fastGridness    = true;       % cached sparse rotation matrices + one-pass annulus correlations (gridness step ~4-5x faster, same scores)
prms.gridPropsV2.commonMask      = 'none';     % 'ell'/'all': compare all rotations on the same bins - tested, no improvement (fewer bins) 

end

% options.smooth (1,:) {mustBeMember(options.smooth,{'boxcar','adaptive'})} = 'boxcar';
%     options.kernel (1,1) {mustBeNumeric} = 5;
%     options.binSizeSpat (1,1) {mustBeNumeric} = 2.5; % in cm^2
%     options.alpha (1,1) {mustBeNumeric} = 200;
%     options.speedFilterFlagRMaps (1,1) {mustBeNumericOrLogical} = false;  % y/n
%     options.speedFilterLimitLow (1,1) {mustBeNumeric}  = 2.5;
%     options.speedFilterLimitHigh (1,1) {mustBeNumeric}  = 400;
%     options.showWaitBar (1,1) {mustBeNumericOrLogical} = false;
%     options.trimNaNs (1,1) {mustBeNumericOrLogical} = false;
%     options.posOnly (1,1) {mustBeNumericOrLogical} = false;
%     options.cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);

 % options.dirSmoothKern (1,1) {mustBeNumeric} = 5;
 %    options.binSizeDir (1,1) {mustBeNumeric} = 6; % in degrees
 %    options.speedFilterFlagDMaps (1,1) {mustBeNumericOrLogical} = false; 
 %    options.speedFilterLimitLow (1,1) {mustBeNumeric}  = 2.5;
 %    options.speedFilterLimitHigh (1,1) {mustBeNumeric}  = 400;
 %    options.showWaitBar (1,1) {mustBeNumericOrLogical} = false;
 %    options.cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);


  % options.minDwellForEdge (1,1) {mustBeNumeric} = 1; % in s
  %   options.durThrCohRun (1,1) {mustBeNumeric} = 2; % In seconds
  %   options.filtSigmaForRunDir (1,1) {mustBeNumeric}= 3;  % Units=sec. Sigma of the Gaussian filter to pre-filter the data before finding CW and CCW runs. Kernel is 2*sigma in length.
  %   options.durThrJump (1,1) {mustBeNumeric} = 2;  % In seconds
  %   options.gradThrForJumpSmooth (1,1) {mustBeNumeric} = 2.5;
  %   options.runDimension {mustBeNumeric} = [];                   % This is the dimension to run (X=1, Y=2)
  %   options.dirTolerance (1,1) {mustBeNumeric} = 70;  % Tolerance for heading direction, relative to arm direction, for calculating run direction (degrees)
  % 
  %   % linear map params
  %   options.showWaitBar (1,1) {mustBeNumericOrLogical} = false;
  %   options.binSizeLinMaps (1,1) {mustBeNumeric} = 2.5; % bin size (cm^2)
  %   options.smoothFlagLinMaps (1,1) {mustBeNumericOrLogical} = true;  % y/n
  %   options.smoothKernelSD (1,1) {mustBeNumeric} = 2;  % SD, in bins
  %   options.speedFilterFlagLMaps (1,1) {mustBeNumericOrLogical} = true;  % y/n
  %   options.speedFilterLimits (1,2) {mustBeNumeric} = [2.5 400]; % See setting above - lock settings for linear and open field together.
  %   options.remTrackEnds (1,1) {mustBeNumericOrLogical} = false;  % Remove this many bins from each end of the linear track. Set to 0 to remove none.
  %   options.normSort (1,1) {mustBeNumericOrLogical} = true;

    % options.minBinProp (1,1) {mustBeNumeric} = 0.005;
    % options.binSizeSpeed (1,1) {mustBeNumeric} = 2; % in cm/s
    % options.maxSpeed (1,1) {mustBeNumeric} = 40;
    % options.smKernelLength (1,1) {mustBeNumeric} = 0.25; % in ms
    % options.normaliseFR (1,1) {mustBeNumericOrLogical} = false;
    % options.confInt (1,1) {mustBeNumeric} = 95;
    % options.showWaitBar (1,1) {mustBeNumericOrLogical} = false;

