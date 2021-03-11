# scaNpix
# General remarks:
•	This is a [Matlab package](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) to load either dacq or neuropixel data into Matlab for inspection in a GUI and further analysis

•	You will probably need to use at least MATLAB 2017, but I coded most of it in 2019, so no guarantees that all will work with older versions than that

•	If you do find a bug or have a request to improve something, please raise an issue rather than email me about it. 

•	Use this code at your own risk; it will hopefully help you get that Nature paper, but equally likely some bug will screw up your analysis 


# Extra toolboxes required:

•	[Circular statistics tool box](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
 
 
# A) Getting Started


## 1. Create the class object:

We first create the object by grabbing some basic parameters (see also section about the parameter space) and then       typically open a UI dialogue to fetch which data to load. This will initialise the object, but not load any data from disk yet

### Syntax: 
```
obj = dacq(prmsMode, uiFlag);
```
   
### Inputs (type):
  
*	_prmsMode_ (string) controls how parameters for loading and pre-processing are collected
  
    * _‘default’_: use default values (default case when no arguments are supplied)
  
    *	_‘ui’_: opens UI dialogue to let user edit values
  
    *	_‘file’_: load params from file (opens UI dialogue as well)
    
* _uiFlag_ (logical) is optional flag (_default=true_) that if false will skip UI dialogues for data selection
  
    *	if _uiFlag=true_ we’ll use a UI dialogue to collect dataset(s) that are to be loaded (however, we won’t load any data from disk just yet).
    
### Examples:
   ```
   obj = dacq;                   % create object and use default parameters
   obj = dacq(‘ui’);             % create object and collect parameters with UI dialogue
   obj = dacq(‘default’, false); % create object using default parameter and skip UI dialogue for set file selection
   ```
## 2. Load the data from disk:
### Syntax:
```
obj.load(loadMode, varargin);
```
### Inputs (type):
   
* _loadMode_ (cell array)
   
     *	controls what part(s) of the data will be loaded into object. Either _{‘all’}_ or any combination of _{‘pos’,’spikes’,’eeg’}_
    
*	_varargin_ (comma separated list of strings)
   
     *	comma separated list of set file names (ommit extensions) to be loaded. Useful when you might want to reload only some particular part of the data
    
### Examples:
```
obj.load;                           % load all trials indicated in obj.trialNames. Types of data to be loaded will be collected with UI dialogue
obj.load({‘pos’,’eeg’});            % load position and eeg data for all trials indicated in obj.trialNames.
obj.load({‘all},’SomeSetFileName’); % load all types of DACQ data for trial ’SomeSetFileName’
```

## 3. Add meta data:
  
Usually we want to add some metadata to be associated with data in object that can’t be directly retrieved from DACQ data files (e.g. age of animal). In the object this will be stored in _obj.metaData_ as a scalar structure with values for each trial (so even global metadata like age will need to be set for each trial). Values are stored in cell arrays to allow for flexibility in input type. 
  
### Syntax:  
```
obj.addMetaData(name, varargin);
```
### Inputs (type):
*	_name_ (string/cell array of strings)
  
     *	Name(s) of field(s) to be added to _obj.metaData_. Note that fields in _obj.defaultMetaDataFields_ are initialised upon loading data into object. 
    
* _varargin_ - optional argument to add values for field(s). 

     * Note that you need to supply values for each trial in object. If a single field is added to obj with just single trial, just supply value directly. For other cases (multiple trials and/or fields) format is _{value field1, value field2,…,value fieldN}_ (multiple trials, but single field) or _{{{values field1},{values field2},…,{ values fieldN}}_ (multiple trials and multiple fields).

### Examples:         
```
obj.addMetaData;                                                         % open UI dialogue to indicate which field and value(s) to add
obj.addMetaData(name);                                                   % initialise metadata field obj.metaData.(name)  with empty cells (e.g. default fields)
obj.addMetaData(name,{values for obj.metaData.(name)-1 for each trial}); % add obj.metaData.(name) to obj and set values=varargin
```

## 4. Make rate maps
  
Making rate maps comes with a fairly large number of possible parameters for which the default values are stored in defaultParamsRateMap.m. If you prefer to use your own combination here is how to (best) do it: 

* Change values as you like in _obj.mapParams.(name)_ and then run _obj.saveParams(‘maps’)_ to store your version as a .mat file in _'PathOnDisk\@dacq\files\'_. 

* Set _obj.params('myRateMapParams')='FilenameYouStoredRateMapParamsAs'_ (no extention). 

* Now when you make maps, we’ll use your version unless specified differently (even though that might sound cumbersome, this makes reverting to defaults straight-forward and you only have to do this once if you also change _obj.params('myRateMapParams')='FilenameOfYourDefaults'_ permanently).

Rate maps will be stored in _obj.maps.(mapsType)_ (2D maps) and _obj.linMaps.(mapType)_ (1D maps), respectively. For details of making different types of firing rate maps, see the respective functions: _makeRateMaps.m, makeDirMaps.m, makeLinRMaps.m_

### Syntax:  
```
obj.makeMaps(mapType, trialInd, varargin );
```
### Inputs:

* _mapType_ (string) 
  
    * _‘rate’_ – make standard 2D rate maps (see _makeRateMaps.m_ for more details). This option will also make spike and position maps. Output is a cell(nCells,1) array of rate/spike/pos maps for each trial (note that if you use _‘boxcar’_ smoothing there is only one position map/trial)
    
    *	_‘dir’_ – make standard directional rate maps (see makeDirMaps.m for more details). Output is a cell(nCells,1) array of directional maps for each trial
    
    *	_‘lin’_ – make rate maps for 1D environments. We will linearise position first and then make the maps (see _makeLinRMaps.m_ and _linearisePosData.m_ for more details). The output format is different to 2D maps – we will get a structure with fields:
    
         * _linRate_: for each trial there will be cell(3,1) array with an nCells-by-position array of firing rates for {1,1} all data or CW {2,1} and CCW {3,1} runs, respectively. 

         * _linPos_: same as _linRate_, but for position only (i.e. pos map)

         * _linRateNormed_: cell(3,2) array with e.g. {1,1} same as _linRate{1}_, but maps normalized by their peak rates and ordered according to their position on the track. {1,2} is index which indicates cell order in {1,1} (so _linRate{1} = linRateNormed{1,1}(linRateNormed{1,2}_)  

*	_trialInd_ (double); index for trial(s) for which you want to make maps (default: all trials)
  
*	_varargin_ - optional argument to control parameter selection (see also section about the parameter space)
  
    *	_'load'_: load parameter file from disk
    
    *	_'ui'_: open UI dialogue for parameter collection
    
    *	_prmsStruct_: structure with parameter fields to be changed from defaults
    
    *	_‘Name’, Value_ pairs: comma separated list of ‘Name’, Value pairs

### Examples:        
```
obj.makeMaps;                              % make maps for all trials with default params and open UI dialogue to choose which type of map to make
obj.makeMaps(‘rate’);                      %  make rate maps for all trials with defaults params
obj.makeMaps(‘lin’, 2);                    % make linearised rate maps with defaults params for trial 2
obj.makeMaps(‘dir’,[], 'ui' );             % make directional rate maps for all trials and open UI dialogue to specify parameter space
obj.makeMaps(‘rate’,1,’ smooth’,’boxcar’); % make rate maps for trial 1, using defaults params, except setting prms.smooth=’boxcar’
```            

## 5. Launch GUI elements

### Syntax:
```
obj.startGUI(GUIType);
```

### Inputs:

*	_GUIType_ (string) 

   *	_‘CluCheck’_ – start GUI to inspect properties of individual clusters after spike sorting. E.g. plot all possible cross-correlograms or directly compare properties of 2 clusters (waveform, cross- and auto-correlograms)
   
   *	_‘eegBrowser’_ – start GUI to browse EEG data (Note: Not implemented yet!)

   *	_‘mapGUI’_ – start GUI to look at rate maps and properties of individual cells in more detail 

### Examples:      
```
obj.startGUI;             % open UI dialogue to indicate which GUI to launch
obj.startGUI(‘clucheck’); % start GUI to inspect quality of spike sorting
```


# B) Parameter space:

There are two different parameter spaces that are used within the MATLAB class.

## 1.	General 
Parameters that are used when loading data and doing some basic pre-processing (e.g. position smoothing). Typically, you will not need to change the majority of these and they are stored in obj.params as a [map container](https://www.mathworks.com/help/matlab/map-containers.html). You can access values by using the name of the individual parameter as the key (e.g. _obj.params(‘posFS’)_ will give you the position sample rate and _obj.params.keys_ will give you a list of all parameters in the container).
The default values are generated with _defaultParamsContainer.m_ and you should leave these as they are, but you can save your own version to a file (_obj.saveParams(‘container’)_ can write current map container in object to disk). You should store your parameter file in _'PathOnYourDisk\Scan_v2\@dacq\files\YourFile.mat'_ 

### Full List:
* _pix/m_: pixel/m – leave empty as will be read from set file
*	_ScalePos2PPM_ – scale position data to this pix/m (_default=400_). This is particularly useful for keeping rate map sizes in proportion, if you recorded data across different environments that have a different size and/or pix/m setting for their tracking 
*	_posMaxSpeed_ – speed > posMaxSpeed will be assumed tracking errors and ignored (set to _NaN_); in m/s (_default=4_)
*	_posSmooth_ – smooth position data over this many ms (_default=400_)
*	_posHead_ – relative position of head to headstage LEDs (_default=0.5_) 
*	_posFs_ – position data sampling rate in Hz; leave empty as will be read from pos file (50Hz)
*	_cutFileType_ – type of cut file, i.e. ‘cut’ (Tint) or ‘clu’ (KlustaKwik); _default=’cut’_
*	_cutTag1_ – cut file tag that follows base filename but precedes \_tetrodeN in filename (_default=’’_); this is something historic (and idiosyncratic) for the data of the original pup replay study, so chances are you can just ignore this 
*	_cutTag2_ – cut file tag that follows \_tetrodeN in filename (_default=’’_)
*	_eegFs_ – low sampling rate EEG data sampling rate in Hz; leave empty as will be read from eeg file (250Hz)
*	_egfFs_ – high sampling rate EEG data sampling rate in Hz; leave empty as will be read from eeg file (4800Hz)
*	_defaultDir_ – default directory for UI dialogues where to look for things, e.g. data (_default='PathToThe@dacqCodeOnYourDisk'_)
*	_myRateMapParams_ – _'FileNameOfYourRateMapParams.mat'_ (_default=’’_)
                
## 2.	Parameters for rate maps. 
These are stored as a scalar MATLAB structure in _obj.mapParams_ (note that _obj.mapParams_ is a hidden property) and the default values are generated with _defaultParamsRateMaps.m_ – again do not edit anything in there. 
You can edit these parameters on the fly when generating different kinds of rate maps.
If you want to use your own custom values by default you should edit them within object and then save them to disk using _obj.saveParams(‘maps’)_ to _PathOnYourDisk\Scan_v2\@dacq\files\YourFile.mat_. Then go and set _obj.params(‘myRateMapParams’)='YourFile’_. 

### Full list:

* General params:
   * _speedFilterLimitLow_ – lower limit for speed in cm/s (_default=2.5_)  
   * _speedFilterLimitHigh_ – upper limit for speed in cm/s (_default=400_)

* 2D rate maps:
   * _speedFilterFlagRMaps_ – logical flag if speed filtering for position data is invoked (_default = true_)
   * _binSizeSpat_ – bin size for spatial rate maps in cm<sup>2</sup> (_default=2.5_) 
   * _smooth_ – type of smoothing; _'adaptive'_ (default) or _'boxcar’_ 
   * _alpha_ – alpha parameter for adaptive smoothing in seconds (_default=200_; usually shouldn’t be changed)
   * _mapSize_ – final size of map (xy) in bins (_default=[ ]_); if empty we will use the max extent of sampled positions, but if you had a badly sampled environment you could reconstruct the sampled space in relation to the actual size (you will need a record of e.g. coordinates of diagonal opposite corners in a rectangular environment)
   
* directional maps:
   * _speedFilterFlagDMaps_ – logical flag if speed filtering for position data is invoked (_default=true_)
   * _binSizeDir_ - bin size for directional maps in degrees (_default=6_)
   * _dirSmoothKern_ – size of smoothing kernel for directional maps in degrees (_default=5_) 
   
* linear rate maps:
   * _trackType_ – (_default=[]_);
   * _binSizeLinMaps_ – bin size for linear rate maps in cm (_default=2.5_)
   * _smoothFlagLinMaps_ – logical flag if maps should be smoothed (_default=true_)
   * _smoothKernelSD_ – SD of Gaussian smoothing kernel in bins (_default=2_). Kernel is 5\*SD in length (should we make this smaller?)
   * _speedFilterFlagLMaps_ – logical flag if speed filtering for position data is invoked (_default=true_)
   * _posIsCircular_ – logical flag position data is assumed to be circular, as e.g. on square track (_default=false_)
   * _remTrackEnds_ – set this many bins to NaN at each end of the linear track (_default=0_). Don’t use for square track data
   * _normSort_ – make normalized and sorted (by rate map peak on track) map array (_default=1_)

* Parameters for linearisation of position data: 
   * _trackLength_ – track length in pixels (_default=[ ]_); as will differ for each type of track). For square track it should be length for 1 arm only. This is crucial to make rate map size match across datasets.
   * _minDwellForEdge_ – minimum dwell of animal in bin at edge of environment in seconds (_default=1_)
   * _durThrCohRun_ – threshold for minimum duration of run in one direction in seconds (_default=2_); set to 0 if you don’t want to remove position data of run periods < threshold
   * _filtSigmaForRunDir_ – SD of the Gaussian filter to pre-filter the data before finding CW and CCW runs in seconds (default=3); Kernel is 2*SD in length.
   * _durThrJump_ – threshold for short periods of change of running direction that will be discounted if gradient < gradThrForJumpSmooth; in seconds (_default=2_)
   * _gradThrForJumpSmooth_ - gradient of the smoothed linear positions in the jump window in cm/s (_default=2.5_);
   * _runDimension_ – will be estimated from data if left empty (_default=[ ]_). Only used for linear track data (somewhat historic parameter as could generally be estimated from data)
   * _dirTolerance_ - Tolerance for heading direction in degrees, relative to arm direction, for calculating run direction on track (_default=70_)


# C) Class properties ( property(format) ):

## General:

* _params_ (containers.Map) – map container with general params
* _dataPath_ (char) – _FullPathToDataOnDisk_
* _dataSetName_ (char) – not used at the moment, but thought as unique identifier for dataset (when and if we can store multiple dacq objects together)
* _trialNames_ (string array) – list of filenames in dataset
* _trial_duration_ (double) – array of trial durations in seconds
* _cellNTetN_ (double) – (nCells,3) array of cell_ID (column 1) and tetrode_ID (column 2) ID for cells. 3rd column contains simple numeric index (1:nCells)
* _metaData_ (struct) – store metadata in fields _metaData.(someString)_    

* _setFileData_ (struct) – non-scalar structure that contains various things from _.set_ file with fields:
   * _tracked_spots_ – n of LEDs
   * _xmin_ – min of camera window X
   * _xmax_ – max of camera window X
   * _ymin_ – min of camera window Y
   * _ymax_ – max of camera window Y
   * _sw_version_ – software version
   * _trial_time_ – start of recording in time of day (as set on machine)
   * _ADC_fullscale_mv_ – scale max for channels at gain=1 in mV (for USB systems should be 1.5V)
   * _lightBearing_ – direction of LEDs in degrees (up to 4 lights)
   * _colactive_ – logic index of active LEDs (probably only relevant if using multi-colour LED tracking in DACQ)
   * _gains_ – nTetrodes x 4 array of channel gains (Note up to 32 tetrodes (128 channels) possible)
   * _fullscale_ – nTetrodes x 4 array of scale max in µV (Note up to 32 tetrodes (128 channels) possible)
   * _eeg_channel_ – nEEGs x 1 array of channels that EEGs were recorded from
   * _eeg_recordingChannel_ – nEEGs x 1 array of channels that were set to EEG in DACQ (this will be same as above if EEG was recorded in mode SIGNAL but differ it was REF)
   * _eeg_slot_ – nEEGs x 1 array of EEG number in DACQ (so .eeg, .eeg2, … , .eegN)
   * _eeg_scalemax_ – nEEGs x 1 array of scale max for EEG channels
   * _eeg_filter_ – nEEGs x 1 array of filter type for EEG (0=DIRECT, 1=DIRECT+NOTCH, 2=HIGHPASS, 3=LOWPASS, 4=LOWPASS+NOTCH)
   * _eeg_filtresp_ – filter type of user defined filter (lowpass, highpass, bandpass, bandstop)
   * _eeg_filtkind_ – filter kind for user defined filter (most likely Butterworth)
   * _eeg_filtfreq1_ – lower bound for user defined bandpass filter (_default=300Hz_)
   * _eeg_filtfreq2_ – upper bound for user defined bandpass filter (_default=7kHz_)
   * _eeg_filtripple_ – mystery parameter of user defined filter (should be left at 0.1 according to manual) 
   
* _posData_ (struct) – scalar structure with position data, with fields:
   * _XY_raw_ – cell arrays of raw animal position in pixels (xy-coordinates)
   * _XY_ – cell arrays of binned (i.e. integers) animal position in pixels (xy-coordinates)
   * _direction_ – cell arrays of head direction of animal in degrees
   * _speed_ – cell arrays of running speed of animal in cm/s
   * _linXY_ – cell arrays of linerised position (will only be generated when making linear rate maps)
   
* spikeData (struct) – scalar structure with spike data, with fields:
   * _spk_times_ – cell arrays of spike times (in seconds)
   * _spk_waveforms_ – cell arrays of waveforms (in µV); format for each cell is nSpikes-by-nSample-by-nChannel (i.e. nSpikesx50x4) 

* _eegData_ (struct) – scalar structure with eeg data, with fields:
   * _rawEEG_ – cell arrays of low sample rate (250Hz) EEG data (in µV)
   * _rawEGF_ – cell arrays of high sample rate (4800Hz) EEG data (in µV)
   * _eegTet_ – cell array of tetrode IDs EEGs were recorded from  
    
* _maps_ (struct) – scalar structure with rate maps
   * _rate_ – cell arrays of standard 2D rate maps
   * _spike_ – cell arrays of 2D spike maps
   * _pos_ – cell arrays of 2D position maps
   * _dir_ – cell arrays of directional maps
   
* _linMaps_ (struct) – scalar structure with linearised rate maps
   * _linRate_ – cell(3,1) arrays of linearised rate maps. Format is full rate maps ({1}), rate map for CW ({2}) and CCW ({3}) runs
   * _linPos_ – cell arrays of linearised position maps. Format is same as for linRate
   * _linRateNormed_ – cell(3,2) arrays of linearised rate maps, normalised by their peak rates and ordered according to their position on the track. Format of first column is as for linRate and in 2nd column you have the index, such that:
_linRate{1}=linRateNormed{1,1}(linRateNormed{1,2}_)

## Hidden properties:

* _uniqueCellIndex_ (logical) – logical index to indicate unique cells (_default=true(nCells,1)_). This is used when calling _obj.findCrossRecCells_ which will set index for those cells which are deemed duplicates to _false_. Then you can either remove them from object or exclude them from further analyses using the index 
* _fields2spare_ (cell array) – cell array of fields that will not be changed when e.g. deleting data from object. Typically, fields that contain only 1 value/dataset (_default={'params','dataPath','dataSetName','cellNTetN'}_)
* _defaultMetaDataFields_ (cell array) – metadata fields that will initialised by default when loading data into object. Will be initialized as empty cells (_default={'age','env'}_)
* _mapParams_ (struct) – stores default rate map params (as in _defaultRateMapParams.m_) 
* _loadFlag_ (logical) – flag to indicate if any data has been loaded into object (_default=false_)
   
 
## D) Class methods:

* to be finished....


