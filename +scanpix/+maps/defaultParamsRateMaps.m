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

% general
prms.speedFilterLimitLow  = 2.5;
prms.speedFilterLimitHigh = 400;
prms.showWaitBar          = false;

% rate maps
prms.speedFilterFlagRMaps = 1;    % y/n
prms.binSizeSpat          = 2.5;  % spatial bin size in cm^2
prms.smooth               = 'adaptive'; % 'boxcar; 'adaptive'
prms.smoothKernel         = 5;    % size smoothing kernel (in bins)
prms.alpha                = 200;  % alpha parameter - don't change
prms.mapSize              = [];

% dir maps
prms.speedFilterFlagDMaps = 1;    % y/n
prms.binSizeDir           = 6;    % in degrees
prms.dirSmoothKern        = 5;    % smooth across that many bins

% lin maps
prms.trackType            = [];
% Parameters for linearisation of position data 
prms.trackLength          = [];   % in pixels; for square track it should be length for 1 arm 
prms.minDwellForEdge      = 1;    % in s
prms.durThrCohRun         = 2;   % In seconds
prms.filtSigmaForRunDir   = 3;    % Units=sec. Sigma of the Gaussian filter to pre-filter the data before finding CW and CCW runs. Kernel is 2*sigma in length.
prms.durThrJump           = 2;    % In seconds
prms.gradThrForJumpSmooth = 2.5;
prms.runDimension         = [];                  % This is the dimension to run (X=1, Y=2); will be estimated if empty (lin track only parameter) 
prms.dirTolerance         = 70;  % Tolerance for heading direction, relative to arm direction, for calculating run direction (degrees, but will switch to radians)
    
% linear map params
prms.binSizeLinMaps       = 2.5; % bin size (cm^2)
prms.smoothFlagLinMaps    = 1;  % y/n
prms.smoothKernelSD       = 2;  % SD, in bins
prms.speedFilterFlagLMaps = 1;  % y/n
prms.posIsCircular        = 0;
prms.remTrackEnds         = 0;  % Remove this many bins from each end of the linear track. Set to 0 to remove none.
prms.normSort             = 1;

% don't change
prms.dirTolerance = prms.dirTolerance * pi/180; % radians

