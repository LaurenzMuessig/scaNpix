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
speedFilterLimitLow              = 2.5;
speedFilterLimitHigh             = 400;
showWaitBar                      = false;

%% rate maps
prms.rate.speedFilterFlagRMaps   = 1;    % y/n
prms.rate.speedFilterLimitLow    = speedFilterLimitLow;
prms.rate.speedFilterLimitHigh   = speedFilterLimitHigh;
prms.rate.binSizeSpat            = 2.5;  % spatial bin size in cm^2
prms.rate.smooth                 = 'adaptive'; % 'boxcar; 'adaptive'
prms.rate.smoothKernel           = 5;    % size smoothing kernel (in bins)
prms.rate.alpha                  = 200;  % alpha parameter - don't change
prms.rate.envSize                = [];
prms.rate.trimNaNs               = false;
prms.rate.posOnly                = false;
prms.rate.posIsFit2Env           = true;
prms.rate.showWaitBar            = showWaitBar;

%% dir maps
prms.dir.speedFilterFlagDMaps    = 1;    % y/n
prms.dir.speedFilterLimitLow     = speedFilterLimitLow;
prms.dir.speedFilterLimitHigh    = speedFilterLimitHigh;
prms.dir.binSizeDir              = 6;    % in degrees
prms.dir.dirSmoothKern           = 5;    % smooth across that many bins
prms.rate.showWaitBar            = showWaitBar;

%% lin maps
prms.lin.trackType               = [];
% Parameters for linearisation of position data 
prms.lin.trackLength             = [];   % in pixels; for square track it should be length for 1 arm 
prms.lin.minDwellForEdge         = 1;    % in s
prms.lin.durThrCohRun            = 2;   % In seconds
prms.lin.filtSigmaForRunDir      = 3;    % Units=sec. Sigma of the Gaussian filter to pre-filter the data before finding CW and CCW runs. Kernel is 2*sigma in length.
prms.lin.durThrJump              = 2;    % In seconds
prms.lin.gradThrForJumpSmooth    = 2.5;
prms.lin.runDimension            = [];                  % This is the dimension to run (X=1, Y=2); will be estimated if empty (lin track only parameter) 
prms.lin.dirTolerance            = 70;  % Tolerance for heading direction, relative to arm direction, for calculating run direction (degrees, but will switch to radians)
    
% linear map params
prms.lin.binSizeLinMaps          = 2.5; % bin size (cm^2)
prms.lin.smoothFlagLinMaps       = 1;  % y/n
prms.lin.smoothKernelSD          = 2;  % SD, in bins
prms.lin.speedFilterFlagLMaps    = 1;  % y/n
prms.lin.speedFilterLimitLow     = speedFilterLimitLow;
prms.lin.speedFilterLimitHigh    = speedFilterLimitHigh;
prms.lin.posIsCircular           = 0;
prms.lin.remTrackEnds            = 0;  % Remove this many bins from each end of the linear track. Set to 0 to remove none.
prms.lin.showWaitBar             = showWaitBar;
% don't change
prms.lin.dirTolerance            = prms.lin.dirTolerance * pi/180; % radians

%% object vector maps
prms.objVect.binSz_dist          = 2.5;   % in cm;  2cm in Høydal et al (2019)
prms.objVect.binSz_dir           = 5;     % in degrees;  5deg in Høydal et al (2019)
prms.objVect.smKernelSz_OV       = 5;
prms.objVect.smSigma_OV          = 2;
prms.objVect.showWaitBar         = showWaitBar;

%% speed maps
prms.speed.minBinProp            = 0.005; % valid speed bins need to contain > prctg of samples of population (0.5%)
prms.speed.speedBinSz            = 2;     % cm/s (2 cm/s)
prms.speed.maxSpeed              = 40;    % cm/s (40 cm/s)
prms.speed.smKernelLength        = 0.25;  % in seconds (250ms)
prms.speed.normaliseFR           = true;  % logical flag
prms.speed.confInt               = 95;    % confidence interval
prms.speed.showWaitBar           = showWaitBar;

%% spatial autocorellograms
prms.sac.method                  = 'moser'; 
prms.sac.removeMinOverlap        = true;
prms.sac.smooth                  = false;
prms.sac.hSize                   = 5;
prms.sac.sigma                   = 1.5;

%% grid props
prms.gridProps.peakDetect        = 'watershed';     % 
prms.gridProps.zScoreThr         = 1;            % 
prms.gridProps.minPeakSz         = 4;    % 
prms.gridProps.plotEllipse       = false;         
prms.gridProps.verbose           = false; 
end

