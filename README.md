# scaNpix
# General remarks:
•	This is a [Matlab package](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) to either load DACQ or neuropixel data into Matlab for inspection in a GUI and further analysis

•	You will need to use at least **MATLAB 2019b** for the code base, but there might be some novel dependencies that I have overlookedintroduced over time that will require a more recent version of Matlab.

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

We first create a data object by grabbing some basic parameters (see also section about the parameter space) and then opening a UI dialogue to fetch which data type (e.g. DACQ or neuropixels) and what files to load. This will initialise the object into the Matlab workspace, but not load any data from disk yet. To select what data to load we first select a parent directory. Then we will get a list of all data files contained within any sub directories and from this list we select which data we want to load (so if your parent directory contains a lot of data files the list will be long). This data will then be treated as one experiment. 

### Syntax: 
```
 obj = scanpix.ephys;
 obj = scanpix.ephys(type);
 obj = scanpix.ephys(type,prmsMode);
 obj = scanpix.ephys(type,prmsMode,setDirFlag);
```
### Inputs (type):
* _type_ (string/character)

     * type of data - _'dacq'_, _'npix'_, _'nexus'_ (not implemented currently) or _'bhave'_ (behavioural data, i.e. no recording data)
  
*	_prmsMode_ (string/character)
   
      *	'default' - uses default parameters defined in _scanpix.helpers.defaultParamsContainer_ (default), 'ui' - opens UI dialogue to select object params; 'file' - load object params from file

 *	_setDirFlag_ (logical)
   
      *	true (default)/false - skip UI set file selection dialogue if false. This stops dialogue pop ups when e.g. batchloading data from disk

Then we can use the class's load method to load the actual data, like so:

```
obj.load(loadMode, varargin);
```
### Inputs (type):
   
* _loadMode_ (cell array)

     * controls what part(s) of the data will be loaded into object. Either _{‘all’}_ or any combination of _{‘pos’,’spikes’,’lfp’}_
     * you can add either of two additional strings: _'nosync'_ to indicate that no sync file with the camera timestamps should be loaded (neuropixel data only; useful when you want to load a concatonated file) or _'reload'_ which will reload the specified data without changing any of the other data (neuropixel data only; this is e.g. useful if you want to load a different version of the spike sorting output that is located inside a different directory) 
  
*	_varargin_ (comma separated list of strings)
   
      *	comma separated list of data file names (ommit extensions) to be loaded. Useful when you might want to reload only some particular part of the data
    
### Examples:
```
obj.load;                           % load all files indicated in obj.trialNames. Types of data to be loaded will be collected with UI dialogue
obj.load({‘pos’,’lfp’});            % load position and eeg data for all files indicated in obj.trialNames.
obj.load({‘all},’SomeDataFileName’); % load all types data for trial ’SomeDataFileName’
```

### Other object methods
• for information on the syntax for using the _obj.methods_, please refer to the descriptions in _scanpix.ephys_
 * _changeParams_: change params of the current object (either object params or map params)
 * _saveParams_: save params of the current object to disk (either object params or map params)
 * _addMetaData_: Add a new field/value pair to the object's metadata
 * _addData_: Add new data to the current object
 * _deleteData_: delete data (trials or cells) from the current object
 * _reorderData_: reorder the sequence of trials within current object
 * _load_: load data
 * _read_histology_: read an xml file with the result of the histology reconstruction (neuropixels only)
 * _deepCopy_: create a deep copy of the current object
 * _addMaps_: add various types of maps (e.g. rate, dir, etc) to the current object
 * _getSpatialProps_: fetch various types of spatial properties (e.g. spatial information, gridness, etc) for the cells in the current object
 * _loadWaves_: load waveforms (neuropixels only as fro dacq these will be loaded during raw data loading)
         


## 3. Do something exciting with the data you loaded into Matlab


### 3A Inspect data in GUI

• The GUI(s) were all made with the [app designer](https://uk.mathworks.com/products/matlab/app-designer.html), so if you want to look at the code or edit something you need to open it in there

• You can start (a) GUI(s) by using a wrapper function (_scanpix.GUI.startGUI_) or by calling it directly (e.g. _mainGUI_)

• When launching the main GUI, we will check if there are any objects matching the data type in the base workspace and ask user if they want to import these, but you can also simply load the data from raw within the GUI or load a GUI state that you saved to disk previously.

• In the GUI you can load and inspect multiple datasets/data from multiple experiments


#### Syntax:
```
scanpix.GUI.startGUI;
scanpix.GUI.startGUI(GUIType);
scanpix.GUI.startGUI(GUIType,dataType);
```
#### Inputs:

*	_GUIType_ (string) 

      *	_‘mainGUI’_ – start main GUI to inspect data
   
      *	_‘lfpBrowser’_ – start GUI to browse EEG data (Note: Not implemented yet!)
 
*	_dataType_ (string) 

      *	_‘dacq’_, _'npix'_ or _'nexus'_

  
or
```
scanpix.GUI.mainGUI;
scanpix.GUI.mainGUI(classType);
```

#### Inputs:

*	_classType_ (string; optional) 

      *	_‘dacq’_, _'npix'_ or _'nexus'_
   

#### Main GUI:

![Picture1](https://user-images.githubusercontent.com/24457903/130802561-6eb1ee3e-7998-4b96-a4a6-0e8986f9db76.png)

#### GUI menu bar:
* File:
     * _Load Data_: Load Datasets from raw, from a GUI state or reload sorting results
       * _from Raw_: load raw data
       * _from GUI state_: load a previous GUI state from disk
       * _batch load_: batch load data (needs formatted spreadsheet, see batch loading above)
       * _reload sorting results_: reload sorting results
     * _Save Data_: Save GUI state; either the currently selected dataset(s) or the full GUI content
       * _Full_: save full GUI content to disk
       * _Current Selection_: only save the GUI with the currently selected datasets
       * _Objects (Current Selection)_: only save the curreently selected datasets as objects (cannot be loaded back into the GUI)
     * _Delete Data_: Delete currently selected dataset(s), trial(s) or cell(s)
     * _Help_: Display help (will display available shortcuts in GUI) 

* Settings:
     * _Set Defaults_
       * _GUI_: Set default parameters for GUI (note that the directory in the GUI defaults is the location on disk where data will be saved to)
       * _General Plotting_: parameters for plots in the _Overview_ tab
       * _Maps_: change params for all map types
       * _Grid cell props calc_: change properties for calculating gridness etc.
       * _Restore GUI defaults_: all aspects of GUI or restore the built-in defaults (this includes map params)
       * _GUI params to objects_: this will push current GUI params (i.e. map params) to all datasets in the GUI
     * _Pause Updates_: Pause updates for any plot(s) on the _Overview_ tab in case updating is slow or plot(s) are irrelevant for user 

* Data:
     * _Add MetaData_: Add meta data to currently selected dataset
       * _add single field_: Add a single new field to currently selected dataset
       * _get nExp_: Add n of exposures for all trial types (currently npix data only). This requires the standard metaData.xml file for each dataset. This will be done for all datasets of an animal that are found in the GUI and will prompt the user to indicate the n of pre-exposures for a given trial type  
     * _Filter_: Filter data according to cluster label (neuropixel data only), minimum number of spikes in any trial or spatial properties. The latter will prompt user to indicate the threshold value, filter direction ('gt'=greater than or 'lt'=lower than) and trial Index. For the trial index either a single numeric index can be used or any number of valid trial indices can be combined by either Matlab AND _&_ or OR _|_ operators (e.g. to filter for cells that pass threshold in the first 3 trials of a given dataset, set trial index to _1&2&3_). Note that currently you cannot combine _&'s_ and _|'s_. Use the _Remove Filter_ option to undo any filtering
     * _Reorder_: reorder datasets in the GUI

* Figures:
     * _DACQ_: nothing here yet
     * _Datasets_: PLot all maps of a certain kind for a given dataset (_PLot All Maps_) or make a huge plot including all map types present in GUI (_Plot All U Got_)
     * _Neuropixels_: Generate data type specific figures (currently only for neuropixel data)
     * _Maps_: Generate time series plots (_Filter by time_ - split trial in time chunks or _Filter by Dir_ - split data by different head directions) or make a plot of grid cell properties (_Plot grid Props_)
     * _Save As PDF_: Save figures generated by GUI as PDFs to disk (either _All_ or _Custom Selection_)


* Analysis:
     * _Spatial Props_: Generate spatial properties for all cells/trials in currently selected data set. Result output is into base work space as a cell array named _datasetName__Property_ with cellID in _output{:,1}_ and values as cell-by-trial array in _output{:,2}_. Note that cells removed by current filter will be ignored for analysis output 

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

#### maps package

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
*	_maxPosInterpolate_ – interpolate over this max distance for between missing positions; in cm (_default=15_)
*	_posHead_ – relative position of head to headstage LEDs (_default=0.5_) 
*	_posFs_ – position data sampling rate in Hz; leave empty as will be read from pos file (50Hz)
*	_cutFileType_ – type of cut file, i.e. ‘cut’ (Tint) or ‘clu’ (KlustaKwik); _default=’cut’_
*	_cutTag1_ – cut file tag that follows base filename but precedes \_tetrodeN in filename (_default=’’_); this is something historic (and idiosyncratic) for the data of the original pup replay study, so chances are you can just ignore this 
*	_cutTag2_ – cut file tag that follows \_tetrodeN in filename (_default=’’_)
*	_APFs_ – data sampling rate in Hz; (48kHz)
*	_loadHighFsLFP_ - true/false; flag to load high sampling rate eeg files
*	_defaultDir_ – default directory for UI dialogues where to look for things, e.g. data (_default='Path/To/The/+scaNpix/Code/On/Your/Disk'_)
*	_myRateMapParams_ – _'FileNameOfYourRateMapParams.mat'_ (_default=’’_)


  
### Full List neuropixel data:

*	_ScalePos2PPM_ – scale position data to this pix/m (_default=400_). This is particularly useful for keeping rate map sizes in proportion, if you recorded data across different environments that have a different size and/or pix/m setting for their tracking 
*	_posMaxSpeed_ – speed > posMaxSpeed will be assumed tracking errors and ignored (set to _NaN_); in m/s (_default=4_)
*	_posSmooth_ – smooth position data over this many s (_default=0.4_)
*	_maxPosInterpolate_ – interpolate over this max distance for between missing positions; in cm (_default=15_)
*	_InterpPos2PosFs_ – true/false; interpolate position data to exact sampling rate (which will be slightly different to exactly 50Hz). This will substantially speed up generating rate maps
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
   * _speedFilterLimits_ – limits for speed filtering in cm/s (_default=[2.5 400]_)  
   * _showWaitBar_ - show waitbar (_default=false_)

* 2D rate maps:
   * _speedFilterFlagRMaps_ – logical flag if speed filtering for position data is invoked (_default = true_)
   * _speedFilterLimitsLow_ – lower limit for speed in cm/s (_default=2.5_); gets set automatically from general params  
   * _speedFilterLimitsHigh_ – upper limit for speed in cm/s (_default=400_); gets set automatically from general params   
   * _binSizeSpat_ – bin size for spatial rate maps in cm<sup>2</sup> (_default=2.5_) 
   * _smooth_ – type of smoothing; _'adaptive'_ (default) or _'boxcar’_
   * _kernel_ – size of boxcar kernel for smoothing (_default=5_)
   * _alpha_ – alpha parameter for adaptive smoothing in seconds (_default=200_; probably shouldn’t be changed)
   * _trimNaNs_ - trim rows or columns in map that are all _NaN_ (_default=false_)
   * _showWaitBar_ - show waitbar (_default=false_); gets set automatically from general params 
   
* directional maps:
   * _speedFilterFlagDMaps_ – logical flag if speed filtering for position data is invoked (_default=true_)
   * _speedFilterLimitsLow_ – lower limit for speed in cm/s (_default=2.5_); gets set automatically from general params  
   * _speedFilterLimitsHigh_ – upper limit for speed in cm/s (_default=400_); gets set automatically from general params  
   * _binSizeDir_ - bin size for directional maps in degrees (_default=6_)
   * _dirSmoothKern_ – size of smoothing kernel for directional maps in degrees (_default=5_)
   * _showWaitBar_ - show waitbar (_default=false_); gets set automatically from general params 
 
   
* linear rate maps:
   * _binSizeLinMaps_ – bin size for linear rate maps in cm (_default=2.5_)
   * _smoothFlagLinMaps_ – logical flag if maps should be smoothed (_default=true_)
   * _smoothKernelSD_ – SD of Gaussian smoothing kernel in bins (_default=2_). Kernel is 5\*SD in length (should we make this smaller?)
   * _speedFilterFlagLMaps_ – logical flag if speed filtering for position data is invoked (_default=true_)
   * _speedFilterLimitsLow_ – lower limit for speed in cm/s (_default=2.5_); gets set automatically from general params  
   * _speedFilterLimitsHigh_ – upper limit for speed in cm/s (_default=400_); gets set automatically from general params 
   * _posIsCircular_ – logical flag position data is assumed to be circular, as e.g. on square track (_default=false_)
   * _remTrackEnds_ – set this many bins to NaN at each end of the linear track (_default=0_). Don’t use for square track data
   * _showWaitBar_ - show waitbar (_default=false_); gets set automatically from general params 
 
* Parameters for linearisation of position data: 
   * _minDwellForEdge_ – minimum dwell of animal in bin at edge of environment in seconds (_default=1_)
   * _durThrCohRun_ – threshold for minimum duration of run in one direction in seconds (_default=2_); set to 0 if you don’t want to remove position data of run periods < threshold
   * _filtSigmaForRunDir_ – SD of the Gaussian filter to pre-filter the data before finding CW and CCW runs in seconds (default=3); Kernel is 2*SD in length.
   * _durThrJump_ – threshold for short periods of change of running direction that will be discounted if gradient < gradThrForJumpSmooth; in seconds (_default=2_)
   * _gradThrForJumpSmooth_ - gradient of the smoothed linear positions in the jump window in cm/s (_default=2.5_);
   * _runDimension_ – will be estimated from data if left empty (_default=[ ]_). Only used for linear track data (somewhat historic parameter as could generally be estimated from data)
   * _dirTolerance_ - Tolerance for heading direction in degrees, relative to arm direction, for calculating run direction on track (_default=70_)

* spatial autocorellograms:
   * _method_ – method to construct autocorellation; '_moser_' or '_barry_' (_default='moser'_)
   * _removeMinOverlap_ – remove bins in sAC where rate map overlap was less than 20 bins (_default=true_)
   * _smooth_ – smooth sAC; not recommended (_default=false_)
   * _hSize_ - size of smoothing kernel (_default=5_)
   * _sigma_ – sigma for smoothing kernel (_default=1.5_)

* grid properties:
   * _binAC_ – bin (or not) sAC before peak detection (_default=true_)
   * _nBinSteps_ – how many bin steps (between -1 and 1) for binning the sAC (_default=21_)
   * _thresh_ – set bins < _thresh_ to NaN (_default=0_)
   * _minPeakSz_ - min number of pixels inside a single peak (_default=8_)
   * _plotEllipse_ – plot the ellipse that was fitted to the inner ring of the sAC (_default=false_)
   * _verbose_ – switch on verbose mode (_default=false_)
   * _legacyMode_ – use the legacy way of computing grid properties (_default=false_)




# C) Class properties ( property(format) ):

## General:

* _params_ (containers.Map) – map container with general params
*  _chanMap_ (struct) – Kilosort channel map (npix data only)
* _dataPath_ (char) – _FullPathToDataOnDisk_ to raw data
* _dataPathSet_ (char) – _FullPathToDataOnDisk_ to sorting output (usually the same as _dataPath_)
* _dataSetName_ (char) – unique identifier for dataset (note that for DACQ this will generate a non-decript name as we don't have a metadata file for Axona data where we could get this information from)
* _trialNames_ (string array) – list of filenames in dataset
* _cellID_ (double) – (nCells,4) array. For DACQ this contains cell_ID (column 1), tetrode_ID (column 2) and simple numeric index (column 3). For neuropixel data this is cluster_ID (column 1), cluster depth (column 2) and channnel closest to centre of mass of cluster (column 3)
* _cellLabel_ (string) – (nCells,1) string array that contains label for clusters from Kilosort ('good' or 'mua') or Phy ('good', 'mua' or 'noise') (npix data only)
* _histo_reconstruct_ (struct) – data for histological reconstruction (npix data only)
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
     * _log_ - log for postion data loading; will give an overview over how good/not so good the tracking data is 
     * _animal_ - animal number
     * _date_ - date of data collection
     * _age_ - age of animal
     * _filename_ - filename for given trial
     * _trialType_ - identifier for trial type (e.g. 'fam' for familiar environment)
     * _duration_ - duration of trial in s
     * _envSize_ size of recording environment in cm
     * _envBorderCoordinates_ - coordinates of recording environment borders in pix.
     * _nLEDs_ - number of LEDs used for tracking
     * _LEDFront_ - which colour LED was at front (if you used 2)
     * _LEDorientation_ - orientation of LEDs with respect to animals' head (if you used 2)
     * _objectPos_ - xy coordinates of objects (if there were any for a given trial)
     * _posFs_ - actual pos sampling rate 
     * _ppm_ - pixel/m – This will contain the final ppm for the position data, i.e. will be different to original when scaling data to common ppm value
     * _ppm__org_ - pixel/m – We store the actual raw data pix/m value here
     * _nChanSort_ - how many channels were used during sorting (usually 383)
     * _nChanTot_ - how many channels in total (usually 385 for SpikeGLX or 384 for OE)
     * _missedSyncPulse_ - list of missing pulses (should ideally be empty)
     * _offset_ - time of first sync pulse in raw data
     * _BonsaiCorruptFlag_ - Bonsai data corrupt; yes/no (can happen in various ways)
     * _PosIsScaled_ - position is scaled to standard ppm; yes/no
   
   
* _posData_ (struct) – scalar structure with position data, with fields:
   * _XY_raw_ – cell arrays of raw LED position data in pixels (xy-coordinates)
   * _XY_ – cell arrays of processed animal position in pixels (xy-coordinates)
   * _direction_ – cell arrays of head direction of animal in degrees
   * _speed_ – cell arrays of running speed of animal in cm/s
   * _linXY_ – cell arrays of linerised position (will only be generated when making linear rate maps)
   * _sampleT_ - sample times of position frames grabbed off Bonsai data. Not really used for anything (npix data only)
   
* spikeData (struct) – scalar structure with spike data, with fields:
   * _spk_times_ – cell arrays of spike times (in seconds)
   * _spk_waveforms_ – cell arrays of waveforms (in µV); format for each cell is nSpikes-by-nSample-by-nChannel (i.e. nSpikesx50x4 for DACQ) 
   * _sampleT_ - timestamps for position frames in neuropixel time. only relevant for npix data

* _lfpData_ (struct) – scalar structure with eeg data, with fields:
   * _lfp_ – cell arrays of low sample rate (250Hz) EEG data for DACQ data or the standard lfp data (2.5kHz)for neuropixel recordings (respectively in µV)
   * _lfpHigh_ – cell arrays of high sample rate (4800Hz) EEG data (in µV) (dacq only) 
   * _lfpTet_ – cell array of tetrode IDs EEGs were recorded from (dacq only) 
    
* _maps_ (struct) – scalar structure with rate maps
   * _rate_ – cell arrays of standard 2D rate maps
   * _spike_ – cell arrays of 2D spike maps
   * _pos_ – cell arrays of 2D position maps
   * _dir_ – cell arrays of directional maps
   * _sACs_ – cell arrays of spatial autocorellograms
   * _lin_ – cell(3,1) arrays of linearised rate maps. Format is full rate maps ({1}), rate map for CW ({2}) and CCW ({3}) runs
   * _linPos_ – cell arrays of linearised position maps. Format is same as for lin
   

## Hidden properties:

* _fileType_ (char) - identifier for data file type (i.e. _.set_ for DACQ and _.ap.bin_ for neuropixel data)
* _fields2spare_ (cell array) – cell array of fields that will not be changed when e.g. deleting data from object. Typically, fields that contain only 1 value/dataset (_default DACQ={'params','dataSetName','cell_ID'}_; _default NPIX={'params','dataSetName','cell_ID','cell_Label'}_)
* _mapParams_ (struct) – stores default rate map params (as in _defaultRateMapParams.m_) 
* _loadFlag_ (logical) – flag to indicate if any data has been loaded into object (_default=false_)
   
 



