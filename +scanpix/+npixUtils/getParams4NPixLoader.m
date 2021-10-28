function prms = getParams4NPixLoader(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Optional inputs/analysis parameters (supply as struct or " ,'fieldname',value, " comma-separated list) :
%
%           prms.loadMode = 'hc'; % 'hc' - hardcoded in prms section; 'ui' - user input
%           % trial names etc - when using 'hc' option
%           prms.rootDir = <directory where data is> - string
%           prms.Tnames = <first and last trial of sequence you want to load> - cell array
%           prms.TNames2exclude = <trials that fall in sequence, but should not be loaded> - cell array
%           % prms for pos files
%           prms.loadPos = 1;
%           prms.ScalePosPPM = 400; %ppm
%           prms.posMaxSpeed = 400; %in cm/s
%           prms.posSmooth = 400; %in ms
%           prms.posHead = 0.5;
%           % prms for tetrode loading
%           prms.loadTet = 1;
%           prms.Tets = 1:8;
%           prms.cutTag2 = 'fin';
%           prms.cutTag1 = '_LTM'; %'_sqTrack'; '';
%           % prms for eeg loading
%           prms.loadEEG = 1;
%           prms.downsampleEGF = 1; %y/n
%           prms.downSampleFact = 4; % 4 is reasonable - will change sr of EGF data to 1200Hz 
%           % prms for map making
%           prms.makeMaps = 1;
%           prms.mapType = 'rate'; % 'dir'; 'rate'
%           %dir maps
%           prms.dirBin = 6; % in degrees
%           prms.dirSmooth = 5;
%           %rate maps
%           prms.placeBin = 10; %cam pix/spatial bin
%           prms.smoothStyle = 'adapt'; % 'adapt' - adaptive smoothing; 'box' - boxcar
%           prms.placeSmooth = 5;  % size smoothing kernel (in bins)
%           prms.adSmooth = 200;  % alpha parameter - don't change
%           %
%           prms.speedFilt = [2.5 400]; % lower/uppper limits cm/s
%           % path scaling (will only work for standard square atm)
%           prms.scalePath = 1; %y/n
%           prms.minOccForEdge = 50; %in pos samples (i.e. 1s)
%           prms.boxExtent = 250; %n of cam pix for 62.5 cm box
%           prms.envCheck = 'famBox'; % for which envs to do path scaling


%% general params
%
prms.paramsFName        = 'tempParams';
%
prms.mode               = 'hc'; % 'hc' - hardcoded in prms section; 'ui' - user input (select trials in dialog)
% trial names etc - when using 'hc' option
prms.dataPath           = ''; % 'Z:\lmuessig\!postDoc\recording_data\CA1\';
prms.Tnames             = {''}; % list of trials to load _famBox _sleepHP _sqTrack
prms.path2RateMapParams = '';

%% params for loading data into class
% prms for pos files
prms.loadPos       = true;
prms.ScalePosPPM   = 400; % scale post data to this ppm (default = 400)
prms.posMaxSpeed   = 4; % in m/s (default = 4 m/s)
prms.posSmooth     = 0.4; % smoothing kernel for pos samples in s (default = 400 ms)
prms.posHead       = 0.5; % head position in relation to LED lights
% prms for tetrode loading
prms.loadSpikes    = true;
% prms for eeg loading
prms.loadEEG       = false;
prms.loadFromPhy   = true;
%% params for path scaling
% path scaling (will only work for standard square atm)
prms.scalePath         = 0; %y/n
prms.triaIndex4Scaling = []; % which trials to do
prms.minOccForEdge     = 1; % in s
prms.boxExtent         = 250; % box size in cam pix if ppm is correct, i.e 250 for a 62.5 cm box @ 400ppm

%% params for maps
% general
prms.speedFilterLimitLow  = 2.5;
prms.speedFilterLimitHigh = 400;
prms.PosFs                = 50;
% rate maps
prms.makeRMaps            = 0;
prms.speedFilterFlagRMaps = 1;    % y/n
prms.binSizeSpat          = 2.5;  % spatial bin size in cm^2
prms.smooth               = 'adaptive'; % 'boxcar; 'adaptive'
prms.smoothKernel         = 5;    % size smoothing kernel (in bins)
prms.alpha                = 200;  % alpha parameter - don't change
prms.mapSize              = [];

% dir maps
prms.makeDirMaps          = 0;
prms.speedFilterFlagDMaps = 1;    % y/n
prms.binSizeDir           = 6;    % in degrees
prms.dirSmoothKern        = 5;    % smooth across that many bins

% lin maps
prms.makeLinMaps          = 0;
% Parameters for linearisation of position data 
prms.minDwellForEdge      = 1;  % in s
prms.durThrCohRun         = 5; % In seconds (10)
prms.filtSigmaForRunDir   = 3;  % Units=sec. Sigma of the Gaussian filter to pre-filter the data before finding CW and CCW runs. Kernel is 2*sigma in length.
prms.durThrJump           = 2;  % In seconds
prms.gradThrForJumpSmooth = 2.5;
prms.runDimension         = [];                  % This is the dimension to run (X=1, Y=2); will be estimated if empty (lin track only parameter) 
prms.dirTolerance         = 70;  % Tolerance for heading direction, relative to arm direction, for calculating run direction.

% linear map params
prms.binSizeLinMaps       = 10; % bin size (pixel) - will give 2.5cm bins
prms.smoothFlagLinMaps    = 1;  % y/n
prms.smoothKernelSD       = 2;  % SD, in bins
prms.speedFilterFlagLMaps = 1;  % y/n
prms.posIsCircular        = 0;
prms.remTrackEnds         = 0;  % Remove this many bins from each end of the linear track. Set to 0 to remove none.
prms.minNSpksInRateMap    = 0;  % min number of spikes for linear rate maps
prms.minPeakRate          = 0;
prms.normSort             = 1;


%% dynamic params
% % dynamic params - these need to be supplied by user - add as many as you
% % like, but make sure they are named valXX (XX needs to go up from 1:1:end)
% prms.metaDataStr   = {'envType','age'};
% prms.val1          = {''};
% prms.val2          = [];


%% overwrite defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prms.dirTolerance         = prms.dirTolerance * pi/180;  % in radians


end

