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
prms.linpos.runDimension          = [];                  % This is the dimension to run (X=1, Y=2); will be estimated if empty (lin track only parameter) 
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
prms.objVect.binSz_dist          = 2.5;   % in cm;  2cm in Høydal et al (2019)
prms.objVect.binSz_dir           = 5;     % in degrees;  5deg in Høydal et al (2019)
prms.objVect.smKernelSz_OV       = 5;
prms.objVect.smSigma_OV          = 2;
prms.objVect.showWaitBar         = showWaitBar;

%% speed maps
prms.speed.minBinProp            = 0.005; % valid speed bins need to contain > prctg of samples of population (0.5%)
prms.speed.binSizeSpeed          = 2;     % cm/s (2 cm/s)
prms.speed.maxSpeed              = 40;    % cm/s (40 cm/s)
prms.speed.smKernelLength        = 0.25;  % in seconds (250ms)
prms.speed.normaliseFR           = true;  % logical flag
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
prms.gridProps.legacyMode        = false; 

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

