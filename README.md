# scaNpix
# General remarks:
•	This is a [Matlab package](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) to either load DACQ or neuropixel data into Matlab for inspection in a GUI and further analysis

•	You will need to use at least **MATLAB 2017b** for the general code base, and **MATLAB 2018b** for the GUI. Note that I coded most of it in MATLAB 2019, so there might be dependencies that I have overlooked.

•	If you do find a bug or have a request to improve something, please [raise an issue in Github](https://docs.github.com/en/github/managing-your-work-on-github/creating-an-issue) rather than email me personally about it. 

•	Use this code at your own risk; it will hopefully help you get that _Nature_ paper, but equally likely some bug will screw up your analysis 


# Extra toolboxes required:

•	[Circular statistics tool box](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
 
 
# A) Getting Started

•	Download/Clone/Photocopy the package and add to your Matlab path (note that Matlab won’t let you add the _+scanpix_ folder itself, but rather only the parent directory)

## 1. Some general remarks about the syntax
•	Functions get called the following way: _scanpix.FunctionName_ or _scanpix.SubPackage.FunctionName_, e.g.if you want to call the function _makeRateMaps_ from the _+maps_ subpackage you would do:
```
scanpix.maps.makeRateMaps(someInput);
```
•	You can use the _Tab_ key to autocomplete for subpackage and/or function name



## 2. Create a class object and load some data:

We first create a data object by grabbing some basic parameters (see also section about the parameter space) and then opening a UI dialogue to fetch which data to load. This will initialise the object into the Matlab workspace, but not load any data from disk yet. To select what data to load we first select a parent directory. Then we will get a list of all data files contained within any sub directories and from this list we select which data we want to load (so if your parent directory contains a lot of data files the list will be long). This data will then be treated as one experiment. 

### Syntax: 
```
obj = scanpix.dacq(someInput);
obj = scanpix.npix(someInput);
```
Then we can use the class's load method to load the actual data, like so:

```
obj.load(loadMode, varargin);
```
### Inputs (type):
   
* _loadMode_ (cell array)

     * controls what part(s) of the data will be loaded into object. Either _{‘all’}_ or any combination of _{‘pos’,’spikes’,’lfp’}_
  
*	_varargin_ (comma separated list of strings)
   
      *	comma separated list of data file names (ommit extensions) to be loaded. Useful when you might want to reload only some particular part of the data
    
### Examples:
```
obj.load;                           % load all files indicated in obj.trialNames. Types of data to be loaded will be collected with UI dialogue
obj.load({‘pos’,’lfp’});            % load position and eeg data for all files indicated in obj.trialNames.
obj.load({‘all},’SomeDataFileName’); % load all types data for trial ’SomeDataFileName’
```


## 3. Do something exciting with the data you loaded into Matlab


### 3A Inspect data in GUI

• The GUI(s) were all made with the [app designer](https://uk.mathworks.com/products/matlab/app-designer.html), so if you want to look at the code or edit something you need to open it in there

• You can start (a) GUI(s) by using a wrapper function (_scanpix.GUI.startGUI_) or by calling it directly (e.g. _mainGUI_)

• When launching the main GUI, we will check if there are any matching data objects in the base workspace and ask user if he wants to import these, but you can also simply load the data from raw within the GUI or load a GUI state that you saved to disk previously.

• In the GUI you can load and inspect multiple datasets/data from multiple experiments


#### Syntax:
```
scanpix.GUI.startGUI(GUIType,classType);
```
or
```
mainGUI(classType);
```

#### Inputs:

*	_GUIType_ (string) 

      *	_‘mainGUI’_ – start main GUI to inspect data
   
      *	_‘lfpBrowser’_ – start GUI to browse EEG data (Note: Not implemented yet!)

*	_classType_ (string) 

      *	_‘dacq’_ – start GUI to inspect DACQ data
   
      *	_‘npix’_ – start GUI to inspect neuropixel data


#### Main GUI:

![Picture1](https://user-images.githubusercontent.com/24457903/130802561-6eb1ee3e-7998-4b96-a4a6-0e8986f9db76.png)

#### GUI menu bar:
* File:
     * _Load Data_ - Load Datasets from raw or a GUI state
     * _Save Data_ - Save GUI state; either the currently selected dataset(s) or the full GUI content
     * _Delete Data_ - Delete currently selected dataset(s), trial(s) or cell(s)
     * _Help_ - Display help (will display available shortcuts in GUI) 

* Settings:
     * _Set Defaults_ - Set default parameters for all aspects of GUI or restore the built-in defaults. Note that the directory in the GUI defaults is the location on disk where data will be saved to   
     * _Pause Updates_ - Pause updates for any plot(s) on the _Overview_ tab in case updating is slow or plot(s) are irrelevant for user 

* Data:
     * _Add MetaData_ - Add meta data to currently selected dataset
     * _Filter_ - Filter data according to cluster label (neuropixel data only), minimum number of spikes in any trial or spatial properties. The latter will prompt user to indicate the threshold value, filter direction ('gt'=greater than or 'lt'=lower than) and trial Index. For the trial index either a single numeric index can be used or any number of valid trial indices can be combined by either Matlab AND _&_ or OR _|_ operators (e.g. to filter for cells that pass threshold in the first 3 trials of a given dataset, set trial index to _1&2&3_). Note that currently you cannot combine _&'s_ and _|'s_. Use the _Remove Filter_ option to undo any filtering  

* Figures:
     * _Save As PDF_ - Save figures generated by GUI as PDFs to disk (either all or custom selection)
     * _Neuropixels or DACQ_ - Generate data type specific figures (currently only for neuropixel data)
     * _Maps_ - Generate time series plots (_Plot time series_) 

* Analysis:
     * _Spatial Props_ - Generate spatial properties for all cells/trials in currently selected data set. Result output is into base work space as a cell array named _datasetName__Property_ with cellID in _output{:,1}_ and values as cell-by-trial array in _output{:,2}_. Note that cells removed by current filter will be ignored for analysis output 

#### GUI shortcuts

* _h_ - display help
* _CTRL+D_ - make Dir Maps
* _CTRL+F_ - save Figures
* _CTRL+L_ - Load data
* _CTRL+M_ - make Rate Maps
* _CTRL+S_ - save GUI state
* _CTRL+F_ - save Figures
* _Up/Down_ - browse through cells
* _Insert/Delete_ - browse through trials
* _PageUp/PageDown_ - browse through datasets
* _Wheelscroll_ - vertical scroll in figures
* _CTRL+Wheelscroll_ - horizontal scroll in figures


#### Batch Loading of data
It's possible to load a batch of datasets into the GUI by supplying paths etc via a spreadsheet. For information on the format of the spreadsheet, see _scanpix.readExpInfo_. 

### 3B Use the code to edit/analyse the data in objects from command window

#### analysis package

• functions/code for data analysis; e.g. calculation of map properties like spatial information etc. as well as cell properties like gridness etc.
```
scanpix.analysis.functionName(someInput);
```

#### dacq/npix Utils maps packages

• functions/code that is specific for data type
```
scanpix.dacqUtils.functionName(someInput);
scanpix.npixUtils.functionName(someInput);
```

#### fxchange package

• functions/code from Matlab file exchange
```
scanpix.fxchange.functionName(someInput);
```

####  GUI packages

• app designer code for GUI(s) as well as a few GUI specific functions. Note that when you debug code in App designer while running the GUI, updates of plots in the GUI tend to become quite sluggish 
```
scanpix.GUI.functionName(someInput);
```

#### helpers package

• functions/code that e.g. helps with data management/processing in objects
```
scanpix.helpers.functionName(someInput);
```

#### maps plot package

• functions/code to generate various types of spatial rate maps
```
scanpix.maps.functionName(someInput);
```

#### plot package

• functions/code to make various types of beautiful plots of your data
```
scanpix.plot.functionName(someInput);
```


# B) Parameter space:

There are two different parameter spaces that are used within _scaNpix_

## 1.	General 
Parameters that are used when loading data and doing some basic pre-processing (e.g. position smoothing). Typically, you will not need to change the majority of these and they are stored in obj.params as a [map container](https://www.mathworks.com/help/matlab/map-containers.html). You can access values by using the name of the individual parameter as the key (e.g. _obj.params(‘posFS’)_ will give you the position sample rate and _obj.params.keys_ will give you a list of all parameters in the container).
The default values are generated with _defaultParamsContainer.m_ and you should leave these as they are, but you can save your own version to a file (_scanpix.helpers.saveParams(obj,‘container’)_ can write the current map container in object to disk). You should store your parameter file in _'PathOnYourDisk\\+scanpix\files\YourFile.mat'_ 

### Full List DACQ:
*	_ScalePos2PPM_ – scale position data to this pix/m (_default=400_). This is particularly useful for keeping rate map sizes in proportion, if you recorded data across different environments that have a different size and/or pix/m setting for their camera setup
*	_posMaxSpeed_ – speed > posMaxSpeed will be assumed tracking errors and ignored (set to _NaN_); in m/s (_default=4_)
*	_posSmooth_ – smooth position data over this many seconds (_default=0.4_)
*	_posHead_ – relative position of head to headstage LEDs (_default=0.5_) 
*	_posFs_ – position data sampling rate in Hz; leave empty as will be read from pos file (50Hz)
*	_cutFileType_ – type of cut file, i.e. ‘cut’ (Tint) or ‘clu’ (KlustaKwik); _default=’cut’_
*	_cutTag1_ – cut file tag that follows base filename but precedes \_tetrodeN in filename (_default=’’_); this is something historic (and idiosyncratic) for the data of the original pup replay study, so chances are you can just ignore this 
*	_cutTag2_ – cut file tag that follows \_tetrodeN in filename (_default=’’_)
*	_lfpFs_ – low sampling rate EEG data sampling rate in Hz; leave empty as will be read from eeg file (250Hz)
*	_lfpHighFs_ – high sampling rate EEG data sampling rate in Hz; leave empty as will be read from eeg file (4800Hz)
*	_defaultDir_ – default directory for UI dialogues where to look for things, e.g. data (_default='Path/To/The/+scaNpix/Code/On/Your/Disk'_)
*	_myRateMapParams_ – _'FileNameOfYourRateMapParams.mat'_ (_default=’’_)
  
  
### Full List neuropixel data:

*	_ScalePos2PPM_ – scale position data to this pix/m (_default=400_). This is particularly useful for keeping rate map sizes in proportion, if you recorded data across different environments that have a different size and/or pix/m setting for their tracking 
*	_posMaxSpeed_ – speed > posMaxSpeed will be assumed tracking errors and ignored (set to _NaN_); in m/s (_default=4_)
*	_posSmooth_ – smooth position data over this many s (_default=0.4_)
*	_posHead_ – relative position of head to headstage LEDs (_default=0.5_) 
*	_posFs_ – position data sampling rate in Hz (_default=50Hz_)
*	_lodFromPhy_ – logical flag to indicate what sorting results to use. If _true_ we'll try Phy otherwise we'll use the raw kilosort results
*	_APFs_ –  sampling rate for single unit neuropixel data (30000Hz)
*	_lfpFs_ – sampling rate for lfp from neuropixel data (2500Hz)
*	_defaultDir_ – default directory for UI dialogues where to look for things, e.g. data (_default='Path/To/The/+scaNpix/Code/On/Your/Disk'_)
*	_myRateMapParams_ – _'FileNameOfYourRateMapParams.mat'_ (_default=’’_)



## 2.	Parameters for rate maps. 
These are stored as a scalar MATLAB structure in _obj.mapParams_ (note that _obj.mapParams_ is a hidden property) and the default values are generated with _scanpix.maps.defaultParamsRateMaps.m_ – again do not edit anything in there. 
You can edit these parameters on the fly when generating different kinds of rate maps.
If you want to use your own custom values by default you should edit them within object and then save them to disk using _scanpix.helpers.saveParams(obj,‘maps’)_ to _PathOnYourDisk/+scanpix/files/YourFile.mat_. Then go and set _obj.params(‘myRateMapParams’)='Path/To/Your/File’_. 

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
   * _binSizeLinMaps_ – bin size for linear rate maps in cm (_default=2.5_)
   * _smoothFlagLinMaps_ – logical flag if maps should be smoothed (_default=true_)
   * _smoothKernelSD_ – SD of Gaussian smoothing kernel in bins (_default=2_). Kernel is 5\*SD in length (should we make this smaller?)
   * _speedFilterFlagLMaps_ – logical flag if speed filtering for position data is invoked (_default=true_)
   * _posIsCircular_ – logical flag position data is assumed to be circular, as e.g. on square track (_default=false_)
   * _remTrackEnds_ – set this many bins to NaN at each end of the linear track (_default=0_). Don’t use for square track data
   * _normSort_ – make normalized and sorted (by rate map peak on track) map array (_default=1_)

* Parameters for linearisation of position data: 
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
*  _chanMap_ (struct) – Kilosort channel map (_scanpix.npix_ only)
* _dataPath_ (char) – _FullPathToDataOnDisk_
* _dataSetName_ (char) – unique identifier for dataset (note that for DACQ this will generate a non-decript name as we don't have a metadata file for Axona data where we could get this information from)
* _trialNames_ (string array) – list of filenames in dataset
* _cellID_ (double) – (nCells,3) array. For DACQ this contains cell_ID (column 1), tetrode_ID (column 2) and simple numeric index (column 3). For neuropixel data this is cluster_ID (column 1), cluster depth (column 2) and channnel closest to centre of mass of cluster (column 3)
* _cellLabel_ (string) – (nCells,1) string array that contains label for clusters from Kilosort ('good' or 'mua') or Phy ('good', 'mua' or 'noise') (_scanpix.npix_ only)
* _metaData_ (struct) – store metadata in fields _metaData.(someString)_    

* _trialMetaData_ (struct) – non-scalar structure that contains various trial specific meta data (from e.g. _.set_ or _.meta_ files :
   * DACQ:  
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
     * _eeg_recordingChannel_ – nEEGs x 1 array of channels that were set to EEG in DACQ (this will be same as above if EEG was recorded in mode SIGNAL but different if it was REF)
     * _eeg_slot_ – nEEGs x 1 array of EEG number in DACQ (so .eeg, .eeg2, … , .eegN)
     * _eeg_scalemax_ – nEEGs x 1 array of scale max for EEG channels
     * _eeg_filter_ – nEEGs x 1 array of filter type for EEG (0=DIRECT, 1=DIRECT+NOTCH, 2=HIGHPASS, 3=LOWPASS, 4=LOWPASS+NOTCH)
     * _eeg_filtresp_ – filter type of user defined filter (lowpass, highpass, bandpass, bandstop)
     * _eeg_filtkind_ – filter kind for user defined filter (most likely Butterworth)
     * _eeg_filtfreq1_ – lower bound for user defined bandpass filter (_default=300Hz_)
     * _eeg_filtfreq2_ – upper bound for user defined bandpass filter (_default=7kHz_)
     * _eeg_filtripple_ – mystery parameter of user defined filter (should be left at 0.1 according to manual) 
     * _ppm_: pixel/m – This will contain the final ppm for the position data, i.e. will be different to original when scaling data to common ppm value
     * _ppm__org_: pixel/m – We store the actual raw data pix/m value here
     * _trackType_ – 'sqtrack or 'linear' (_default=''_);      
     * _trackLength_ – track length in cm (_default=[ ]_); as will differ for each type of track). For square track it should be length for 1 arm only. This is crucial to make rate map size match across datasets.


   
   * Neuropixel data: 
     * _PARAMETERS GO HERE_
     * _ppm_: pixel/m – This will contain the final ppm for the position data, i.e. will be different to original when scaling data to common ppm value
     * _ppm__org_: pixel/m – We store the actual raw data pix/m value here
     * _trackType_ – 'sqtrack or 'linear' (_default=''_);      
     * _trackLength_ – track length in cm (_default=[ ]_); as will differ for each type of track). For square track it should be length for 1 arm only. This is crucial to make rate map size match across datasets.
   
* _posData_ (struct) – scalar structure with position data, with fields:
   * _XY_raw_ – cell arrays of raw animal position in pixels (xy-coordinates)
   * _XY_ – cell arrays of binned (i.e. integers) animal position in pixels (xy-coordinates)
   * _direction_ – cell arrays of head direction of animal in degrees
   * _speed_ – cell arrays of running speed of animal in cm/s
   * _linXY_ – cell arrays of linerised position (will only be generated when making linear rate maps)
   * _sampleT_ - sample times of position frames grabbed off Bonsai data. Not really used for anything (_scanpix.npix_ only)
   
* spikeData (struct) – scalar structure with spike data, with fields:
   * _spk_times_ – cell arrays of spike times (in seconds)
   * _spk_waveforms_ – cell arrays of waveforms (in µV); format for each cell is nSpikes-by-nSample-by-nChannel (i.e. nSpikesx50x4 for DACQ) 
   * _sampleT_ - timestamps for position frames in neuropixel time. only relevant for _scanpix.npix_

* _lfpData_ (struct) – scalar structure with eeg data, with fields:
   * _lfp_ – cell arrays of low sample rate (250Hz) EEG data for DACQ data or the standard lfp data (2.5kHz)for neuropixel recordings (respectively in µV)
   * _lfpHigh_ – cell arrays of high sample rate (4800Hz) EEG data (in µV) (_scanpix.dacq_ only) 
   * _lfpTet_ – cell array of tetrode IDs EEGs were recorded from (_scanpix.dacq_ only) 
    
* _maps_ (struct) – scalar structure with rate maps
   * _rate_ – cell arrays of standard 2D rate maps
   * _spike_ – cell arrays of 2D spike maps
   * _pos_ – cell arrays of 2D position maps
   * _dir_ – cell arrays of directional maps
   * _sACs_ – cell arrays of spatial autocorellograms
   
* _linMaps_ (struct) – scalar structure with linearised rate maps
   * _linRate_ – cell(3,1) arrays of linearised rate maps. Format is full rate maps ({1}), rate map for CW ({2}) and CCW ({3}) runs
   * _linPos_ – cell arrays of linearised position maps. Format is same as for linRate
   * _linRateNormed_ – cell(3,2) arrays of linearised rate maps, normalised by their peak rates and ordered according to their position on the track. Format of first column is as for linRate and in 2nd column you have the index, such that:
_linRate{1}=linRateNormed{1,1}(linRateNormed{1,2}_)

## Hidden properties:

* _fileType_ (char) - identifier for data file type (i.e. _.set_ for DACQ and _.ap.bin_ for neuropixel data)
* _fields2spare_ (cell array) – cell array of fields that will not be changed when e.g. deleting data from object. Typically, fields that contain only 1 value/dataset (_default DACQ={'params','dataSetName','cell_ID'}_; _default NPIX={'params','dataSetName','cell_ID','cell_Label'}_)
* _mapParams_ (struct) – stores default rate map params (as in _defaultRateMapParams.m_) 
* _loadFlag_ (logical) – flag to indicate if any data has been loaded into object (_default=false_)
   
 



