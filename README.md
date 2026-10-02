# scaNpix

## General remarks

* This is a [Matlab package](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) to load **DACQ (Axona) tetrode data**, **Neuropixels (1.0) data** or **purely behavioural (tracking-only) data** into Matlab. You can then look at the data in a GUI or analyse it from the command window.
* All data live in one class (`scanpix.ephys`). The class works the same way for every data type, so analysis code written for one recording system usually runs on the others too.
* **Matlab version:** you need **at least R2020b**, because the code makes heavy use of `arguments` blocks with validators like `mustBeA`/`mustBeText`. Some plotting functions (e.g. `scanpix.plot.plotGridProps`) use `clim`, which needs **R2022a or newer**, so use R2022a+ to be safe. Other newer dependencies may have slipped in that I've overlooked.
* If you find a bug or want something improved, please [raise an issue on GitHub](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue) rather than emailing me.
* Use this code at your own risk. It will hopefully help you get that _Nature_ paper, but some bug could just as easily screw up your analysis.

## Requirements

### Matlab toolboxes
* **Image Processing Toolbox**: rate map smoothing, field detection, grid properties, border score (`imfilter`, `fspecial`, `watershed`, `regionprops`, `bwlabel`, ...)
* **Signal Processing Toolbox**: LFP filtering/theta phase and Neuropixels waveform extraction (`filtfilt`, `butter`, `hilbert`)
* **Statistics and Machine Learning Toolbox**: e.g. `prctile`, `pdist` in shuffling and field detection

### Third party code
* [Circular statistics toolbox](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics): used by e.g. `scanpix.analysis.gridprops` and `scanpix.maps.linearisePosData`
* [npy-matlab](https://github.com/kwikteam/npy-matlab) (`readNPY`): needed to load Kilosort/Phy output for Neuropixels data
* A few file exchange functions ship with the package in `+scanpix/+fxchange` (colour maps, xml reading/writing, `shadedErrorBar`, etc.), so you don't need to install those.

### Optional
* **Python with `scipy` + `matplotlib`**: only needed to export figures as real vector PDFs through the matplotlib route (`scanpix.helpers.saveFigAsPDF(..., 'contentType', 'python')`, see [Saving figures](#saving-figures-as-pdf)).

---

# A) Getting started

* Download/Clone/Photocopy the package and add it to your Matlab path. Matlab won't let you add the `+scanpix` folder itself, only its parent directory.
* `+scanpix/files/` holds the parameter files (`.mat`). Your own parameter files go here too (see [Parameter space](#b-parameter-space)).
* `useful files/` holds a template for the metadata xml file (`masterXML.xml`) that Neuropixels/behavioural data need, plus a GUI manual (`GUI_manual.pptx`).

## 1. General remarks about the syntax

* Call functions as `scanpix.FunctionName` or `scanpix.SubPackage.FunctionName`. For example, to call `makeRateMaps` from the `+maps` subpackage:

```matlab
scanpix.maps.makeRateMaps(someInput);
```

* You can use the _Tab_ key to autocomplete subpackage and function names.
* Most functions have a help header (`help scanpix.maps.makeRateMaps`). Many use Matlab's name-value syntax for optional inputs.

## 2. Expected data layout on disk

### DACQ
* The standard Axona files for a trial (`.set`, `.pos`, `.N` tetrode files, `.eeg`/`.eegN`, `.egf`/`.egfN`) and the cut files (`.cut` from Tint or `.clu` from KlustaKwik; see `cutFileType`, `cutTag1`, `cutTag2` [below](#full-list-dacq)).
* **Optional:** an xml file named `meta*_<trialName>.xml` in the same folder. Every field in it (except comments) is added to `obj.trialMetaData` for that trial. This is the way to set e.g. `envSize` and `trialType` for DACQ data, and the code needs these to fit positions to the environment (see `scalePos2CamWin`).

### Neuropixels
Each trial lives in its own folder, which contains:
* the raw data (`*.ap.bin` + `*.ap.meta`, and `*.lf.bin` + `*.lf.meta` if you want LFP)
* a metadata xml file (`*meta*.xml`) with trial info (animal, date, age, trial type, duration, environment size and borders, LED setup, object positions, ...). See `useful files/masterXML.xml` for the template.
* a Kilosort channel map (`*kiloSortChanMap*`)
* the tracking data from Bonsai as a `.csv` file in a subfolder `trackingData/`
* the sync TTLs (`*TTL.txt` or `syncTTLs.mat`). If neither exists, the TTLs are extracted from the sync channel of the raw data.
* the spike sorting output from Kilosort (`spike_times.npy`, `spike_clusters.npy`, `templates.npy`, `cluster_KSLabel.tsv`, ...) and, if curated, from Phy (`cluster_info.tsv`). If the files aren't in the data folder, you're asked where they are; that location is stored in `obj.dataPathSort`.
* **Optional:** a `*hist*.xml` file with the histological reconstruction of the probe track (see `obj.read_histology`)

### Behavioural data (`'bhave'`)
* A Bonsai tracking `.csv` file and an `.xml` metadata file per trial. No ephys data.

## 3. Create a class object and load some data

First we create a data object. This grabs some basic parameters (see the [parameter space](#b-parameter-space) section) and opens a UI dialogue that asks which data type to use (e.g. DACQ or Neuropixels) and which files to load. The object now exists in the Matlab workspace, but no data has been loaded from disk yet.

To pick the data, you first select a parent directory. You then get a list of all data files in any of its subdirectories and select the ones you want to load. If your parent directory holds many data files, this list will be long. All the data you select are treated as one experiment (dataset).

### Syntax

```matlab
obj = scanpix.ephys;
obj = scanpix.ephys(type);
obj = scanpix.ephys(type, prmsMode);
obj = scanpix.ephys(type, prmsMode, setDirFlag);
```

### Inputs
* _type_ (char)
  * type of data: `'dacq'`, `'npix'`, `'bhave'` (behavioural data, i.e. no recording data) or `'nexus'` (not implemented currently)
* _prmsMode_ (char)
  * `'default'` uses the default parameters from `scanpix.helpers.defaultParamsContainer` (default)
  * `'ui'` opens a UI dialogue to set the object params
  * `'file'` loads the object params from a file
* _setDirFlag_ (logical)
  * true (default) / false. If false, the file selection dialogue is skipped. This stops dialogues from popping up when you e.g. batch load data from disk.

Then use the class's `load` method to load the actual data:

```matlab
obj.load(loadMode, varargin);
```

### Inputs
* _loadMode_ (cell array)
  * controls which part(s) of the data are loaded into the object: either `{'all'}` (= position + spikes) or any combination of `{'pos','spikes','lfp'}`
  * you can add either of two extra strings:
    * `'nosync'`: don't load a sync file with the camera timestamps (Neuropixels only; useful when you load a concatenated file)
    * `'reload'`: reload the specified data without changing any of the other data (Neuropixels only; useful e.g. for loading a different version of the spike sorting output from another directory)
* _varargin_ (comma separated list of strings)
  * names of the data files to load (omit extensions). Useful when you only want to (re)load a particular part of the data.

Note that `'all'` does **not** load LFP. Ask for `'lfp'` explicitly. For Neuropixels data, `obj.lfpParams` sets which channels are loaded (see [LFP parameters](#3-parameters-for-loading-neuropixels-lfp)).

### Examples

```matlab
obj.load;                              % load all files listed in obj.trialNames; a UI dialogue asks which data types to load
obj.load({'pos','lfp'});               % load position and LFP data for all files in obj.trialNames
obj.load({'all'},'SomeDataFileName');  % load position + spike data for trial 'SomeDataFileName'
```

## 4. Loading data programmatically (no UI)

### Single datasets: `scanpix.objLoader`

Creates an object and loads data from raw in one go, with no dialogues. Only position and spike data are loaded by default.

```matlab
obj = scanpix.objLoader(objType, dataPath);
obj = scanpix.objLoader(objType, dataPath, metaData);
obj = scanpix.objLoader(objType, dataPath, metaData, 'paramsobj', prmsMap, 'paramsmap', mapPrms, 'loadlfp', true);
```

* _objType_: `'npix'`, `'dacq'` or `'bhave'`
* _dataPath_: cell array of full paths to the data files (e.g. `{'X:\data\r123\trial1.ap.bin', ...}`). If you pass a concatenated `.dat` file, the object is flagged as concatenated (`obj.isConcat`).
* _metaData_ (optional): nx2 cell array of `{fieldName, values}` that is added to `obj.trialMetaData`. For DACQ, `cutTag1`/`cutTag2` can be passed here as well.
* name-value pairs:
  * `paramsobj`: containers.Map with the object params (default: `scanpix.helpers.defaultParamsContainer(objType)`)
  * `paramsmap`: struct with the map params (default: `scanpix.maps.defaultParamsRateMaps`)
  * `paramslfp`: struct with the Neuropixels LFP loading params (default: `scanpix.helpers.defaultParamsLFP`)
  * `loadpos` / `loadspikes` / `loadlfp`: logical flags (defaults: true / true / false)

### Batch loading: `scanpix.batchLoader`

Loads a whole batch of datasets based on a spreadsheet (a wrapper around `scanpix.objLoader`):

```matlab
objData = scanpix.batchLoader(cribSheetPath, method, objType, addMeta, Name, Value);
```

* _cribSheetPath_: path to an Excel sheet. It must have the columns `filePath`, `animal`, `trialName` and `experiment_ID`. You can add any other columns, e.g. age, and pull them into the object metadata with `addMeta`. See `scanpix.helpers.readExpInfo` for details on the format.
* _method_: `'exp'` (group rows into experiments by animal + `experiment_ID`; default), `'single'` (every row is its own dataset), `'singleExp'` or `'singleTrial'`
* _objType_: `'npix'` (default), `'dacq'` or `'bhave'`
* _addMeta_: cell array of column names from the sheet to add as metadata
* name-value pairs: same as for `scanpix.objLoader`

The output is a cell array of objects. Datasets that fail to load are reported, and loading carries on with the rest.

Other helpers that are handy for building analysis loops:
* `scanpix.helpers.getTrialSequence` / `batchGetTrialSequence`: find the trials of given types in a raw data directory (from the xml metadata)
* `scanpix.helpers.matchTrialSeq2Pattern`: match a pattern of trial types to the trials run in a dataset

## 5. Object methods

For details on the syntax, see the help for each method in `scanpix.ephys` (e.g. `help scanpix.ephys.addMaps`).

| Method | Description |
|---|---|
| `load` | load data (see above) |
| `changeParams` | change the params of the current object (object params or map params) |
| `saveParams` | save the params of the current object to disk (`'container'` = object params, `'maps'` = map params) |
| `addMetaData` | add a new field/value pair to `obj.trialMetaData` |
| `addData` | add more trials to the current object (existing maps and metadata fields are removed and need to be regenerated) |
| `deleteData` | delete trials or cells from the current object |
| `reorderData` | reorder the trials in the current object |
| `truncateData` | cut a trial down to a time window `[tStart tEnd]` (e.g. when the headstage came unplugged halfway through). Spike times, position and LFP are re-referenced so the trial starts at t=0 again. |
| `read_histology` | read an xml file with the result of the histological reconstruction (Neuropixels only; called automatically during `load`) |
| `deepCopy` | create a deep copy of the current object (the class is a `handle`, so `obj2 = obj` does **not** copy!) |
| `addMaps` | add maps to the object: `'pos'`, `'rate'`, `'dir'`, `'lin'`, `'sac'` (spatial autocorrelograms), `'objVect'` (object vector maps) or `'speed'` |
| `getSpatialProps` | compute spatial properties for all cells: `'si'` (spatial info), `'rv'` (Rayleigh vector length), `'gridness'` (gridness, wavelength, orientation, both regular and from the ellipse fit), `'bs'` (border score) |
| `loadWaves` | load waveforms (Neuropixels only; DACQ waveforms are loaded together with the spikes) |

### Example workflow

```matlab
obj = scanpix.ephys('npix','default');      % select data in UI
obj.load({'pos','spikes'});
obj.addMaps('rate');                          % rate maps for all trials
obj.addMaps('sac');                           % needs rate maps first
gridProps = obj.getSpatialProps('gridness');  % nCells x 6 x nTrials
SI        = obj.getSpatialProps('si');        % nCells x nTrials
scanpix.plot.plotRateMap(obj.maps.rate{1}{1});
```

## 6. Do something exciting with the data you loaded into Matlab

### 6A Inspect data in the GUI

* The GUIs were all made with [App Designer](https://uk.mathworks.com/products/matlab/app-designer.html). If you want to look at or edit the code, you need to open them there.
* You can start a GUI with the wrapper function `scanpix.GUI.startGUI` or by calling it directly (e.g. `scanpix.GUI.mainGUI`).
* When the main GUI launches, it checks the base workspace for objects of the matching data type and asks if you want to import them. You can also load data from raw within the GUI, or load a GUI state you saved to disk earlier.
* You can load and inspect several datasets (from different experiments) in the GUI at the same time.

#### Syntax

```matlab
scanpix.GUI.startGUI;
scanpix.GUI.startGUI(GUIType);
scanpix.GUI.startGUI(GUIType, dataType);
```

#### Inputs
* _GUIType_ (char)
  * `'main'`: start the main GUI to inspect data (default)
  * `'lfpBrowser'`: GUI to browse LFP data (**not implemented yet!**)
* _dataType_ (char)
  * `'npix'` (default) or `'dacq'`

or

```matlab
scanpix.GUI.mainGUI;
scanpix.GUI.mainGUI(classType);  % classType: 'dacq' or 'npix' (optional)
```

#### Main GUI

![Picture1](https://user-images.githubusercontent.com/24457903/130802561-6eb1ee3e-7998-4b96-a4a6-0e8986f9db76.png)

#### GUI tabs
* **Overview**: browse datasets/trials/cells and look at rate maps, directional maps and other map types (spatial ACs, object vector maps, speed maps), waveforms etc. Make maps for the current dataset with the _Make Maps_ / _Make Plots_ buttons.
* **Compare Cells**: plot temporal cross-correlograms (_Plot CG Selection_), waveforms (_Plot WF Selection_) and the cluster space (_Plot Cluster Space_) for a selection of cells
* **Linearise It**: linearise position data and make linear rate maps for linear/square track data
* **VTC analysis**: vector trace cell analysis (baseline/probe/post-probe trials, vector maps, trace and overlap scores, batch scoring and export; see the `+vtc` package)
* **LFP browser**: placeholder, not implemented yet

#### GUI menu bar
* **File**
  * _Load Data_: load datasets
    * _from Raw_: load raw data
    * _from GUI state_: load a previous GUI state from disk
    * _from Objects_: import a cell array of `scanpix.ephys` objects from an `objData_*.mat` file (e.g. one saved with _Save Data > Objects (Current Selection)_)
    * _batch load_: batch load data (needs a formatted spreadsheet, see [batch loading](#batch-loading-scanpixbatchloader))
    * _reload sorting results_: reload the spike sorting results (e.g. after more curation in Phy)
  * _Save Data_: save the GUI state, either the currently selected dataset(s) or the full GUI content
    * _Full_: save the full GUI content to disk
    * _Current Selection_: only save the GUI with the currently selected datasets
    * _Objects (Current Selection)_: only save the currently selected datasets as objects (these can't be loaded back as a GUI state, but you can import them with _Load Data > from Objects_)
  * _Delete Data_: delete the currently selected _Dataset(s)_, _Trial(s)_ or _Cell(s)_
  * _Help_: show help (the available shortcuts in the GUI)
* **Settings**
  * _Set Defaults_
    * _GUI_: set the default parameters for the GUI. The directory in the GUI defaults is where data are saved to.
    * _General Plotting_: parameters for the plots in the _Overview_ tab
    * _Maps_: change the params for each map type (_Rate_, _Dir_, _Linear_, _ObjectVect_, _Speed_, _sAC_)
    * _Grid cell props calc._: change the properties for calculating gridness etc.
    * _Restore GUI defaults_: restore the built-in defaults for all aspects of the GUI (including the map params)
    * _GUI params to objects_: push the current GUI params (i.e. map params) to all datasets in the GUI
  * _Pause Updates_ (_Freeze_ / _UnFreeze_): pause updates for any plot(s) on the _Overview_ tab, in case updating is slow or the plot(s) don't matter to you
* **Data**
  * _Add MetaData_: add metadata to the currently selected dataset
    * _add single field_: add a single new field to the currently selected dataset
    * _get nExp_: add the number of exposures for all trial types (Neuropixels only). This needs the standard metadata xml file for each dataset. It runs for all datasets of an animal found in the GUI and asks you for the number of pre-exposures for a given trial type.
  * _Edit dataset name_: rename the current dataset
  * _Filter_: filter the cells in a dataset
    * _By Clust Label_: by cluster label (Neuropixels only)
    * _By N Spikes_: by a minimum number of spikes in any trial
    * _By Spatial Props_: by _spatial Info_, _gridness_ or _rayleigh vector_. You're asked for the threshold value, the filter direction (`'gt'` = greater than or `'lt'` = lower than) and the trial index. The trial index is either a single number, or several valid trial indices combined with Matlab AND `&` or OR `|` (e.g. to keep cells that pass the threshold in the first 3 trials of a dataset, set the trial index to `1&2&3`). You currently can't mix `&` and `|`.
    * _By Selection_: keep a manual selection of cells
    * _By Shuffle_: by thresholds from shuffled data (see the `+shuffle` package)
    * _Custom_: custom filter
    * _Remove Filter_: undo any filtering
  * _Reorder_: reorder the _Datasets_ in the GUI
  * _Truncate_: truncate a trial to a time window (see `obj.truncateData`)
* **Figures**
  * _DACQ_: nothing here yet
  * _Dataset_: plot all maps of a certain kind for a dataset (_Plot All Maps_), or make one huge plot with all map types in the GUI (_Plot All U Got_; better go and get a coffee...)
  * _Neuropixels_: _plot Waveform_ and _plot Wave Dist on Probe_ (distribution of the units' waveforms along the shank)
  * _Maps_: time/direction/speed series plots (_Filter by time_ splits the trial into time chunks, _Filter by Dir_ splits the data by head direction, _Filter by speed_ splits it by running speed), or a plot of grid cell properties (_Plot Grid Props_)
  * _Save as PDF_: save figures generated by the GUI to disk as PDFs (_All_ or _Custom Selection_)
* **Analysis**
  * _Spatial props_: compute spatial properties (_grid properties_, _spatial info_, _rayleigh vector length_, _border score_, _spatial correlation_) for all cells/trials in the currently selected dataset. The result goes to the base workspace as a cell array named `datasetName_Property`, with the cellID in `output{:,1}` and the values as a cell-by-trial array in `output{:,2}`. Cells removed by the current filter are left out.

#### GUI shortcuts
* _h_: display help
* _CTRL+D_: make Dir Maps
* _CTRL+R_: make Rate Maps
* _CTRL+F_: save Figures
* _CTRL+L_: Load data
* _CTRL+S_: save GUI state
* _Up/Down_: browse through cells
* _Insert/Delete_: browse through trials
* _PageUp/PageDown_: browse through datasets
* _Wheelscroll_: vertical scroll in figures
* _CTRL+Wheelscroll_: horizontal scroll in figures

### 6B Use the code to edit/analyse the data in objects from the command window

#### analysis package
Functions for data analysis: properties of maps (spatial information, sparsity, gridness, border score, spatial correlations, field detection, ...) and of cells (waveform properties, temporal auto/cross-correlograms, STTC, speed tuning, ...).

```matlab
scanpix.analysis.functionName(someInput);
```

Some highlights:
* `gridprops`: gridness, wavelength, orientation and offset from a spatial AC (with an optional ellipse fit), `selectBestGridProps`, `computeIURatio`, `getPhaseOffset`
* `sortModules`: sort grid cells into modules
* `fieldDetect`: find fields/peaks in any rate map or spatial AC
* `spatial_info`, `getSparsity`, `getBorderScore`, `rayleighVect`, `getMeanRate`
* `spatialCorrelation` (fully vectorised), `spatialCrosscorr`
* `spk_crosscorr`, `get1stMomentAC`, `getSTTC` (spike time tiling coefficient), `getWaveFormProps`
* `classifySpeedTuning`, `getCIProportions`
* `obj2table`: convert a `scanpix.ephys` object into a table (one row per cell; can add maps, grid props, waveform props and LFP). Handy for pooling data across datasets.

#### dacqUtils / npixUtils / bhaveUtils packages
Functions specific to a data type (loading raw data, sync, channel maps, waveform extraction, ...).

```matlab
scanpix.dacqUtils.functionName(someInput);
scanpix.npixUtils.functionName(someInput);
scanpix.bhaveUtils.functionName(someInput);
```

#### lfp package
Functions for working with LFP data:
* `lfpFilter`: band pass filter around a peak frequency (zero phase FIR, as in Scan)
* `lfpPowerSpec`: power spectrum; peak frequency, power and signal-to-noise ratio in a band
* `getThetaPhase`: phase per sample and cycle segmentation (removes low power cycles, cycles of the wrong length and phase slips)
* `speedFilter2lfp`: convert a speed filter in position samples into one in LFP samples
* `lfp2uV`: convert Neuropixels LFP (stored as raw `int16` in the object) into µV

```matlab
[lfpUV, t, chans] = scanpix.lfp.lfp2uV(obj, trialInd);
thetaFilt         = scanpix.lfp.lfpFilter(lfpUV(1,:), 8, obj.trialMetaData(trialInd).lfpFs);
```

#### maps package
Functions to make various types of maps: spatial rate maps, directional, linear, object vector and speed maps. Also time/direction series of maps (`makeMapTimeSeries`, `makeMapDirSeries`), position scaling and linearisation, adaptive smoothing, colour maps.

```matlab
scanpix.maps.functionName(someInput);
```

#### plot package
Functions to make various types of beautiful plots of your data (rate maps, dir maps, linear maps, speed maps, spikes on path, waveforms, cross-correlograms, grid props, scrollable multi-panel plots, ...).

```matlab
scanpix.plot.functionName(someInput);
```

#### shuffle package
Functions to build shuffled distributions (by shifting spike times) and get significance thresholds for cell classification, either per cell or per population (e.g. by age bin).

```matlab
ResShuf = scanpix.shuffle.generateShuffData(objData, 'cell');
thresh  = scanpix.shuffle.getShufThresh(shufVals);
ResT    = scanpix.shuffle.addShufData2Table(ResT, ResShuf, scores);
```

#### vtc package
Functions for analysing vector trace cells (VTCs): vector maps relative to walls/barriers (`makeVectMapBVC`), defining the main field, trace and overlap scores, plotting. The _VTC analysis_ tab in the main GUI uses these functions.

```matlab
scanpix.vtc.functionName(someInput);
```

#### GUI package
App Designer code for the GUI(s), plus a few GUI-specific functions. When you debug code in App Designer while the GUI is running, plot updates in the GUI tend to become quite sluggish.

```matlab
scanpix.GUI.functionName(someInput);
```

#### helpers package
Functions that help with e.g. data management/processing in objects, parameter defaults, reading crib sheets, tables, and saving figures.

```matlab
scanpix.helpers.functionName(someInput);
```

#### fxchange package
Functions from the Matlab file exchange.

```matlab
scanpix.fxchange.functionName(someInput);
```

### Saving figures as PDF

```matlab
scanpix.helpers.saveFigAsPDF(figHandle, 'filename', 'myFig', 'dir', 'C:\figs', 'contentType', 'image');
```

* `'image'` (default): rasterised PDF at `'resolution'` dpi (default 300)
* `'vector'`: Matlab's own vector export. **Caution:** a current MathWorks bug corrupts rasterised content (rate maps etc.) in vector PDF/EPS export.
* `'python'`: rebuilds the figure in matplotlib (`scanpix.helpers.exportViaPython`) to get a real vector PDF with editable text. This needs a Python with `scipy` + `matplotlib`. Set `'pythonExe'` to point to it. Images, lines, areas and text are supported. Polar axes, surfaces, patches, legends, colorbars, scatter and bar plots are skipped with a warning.

Scrollable plots (e.g. from `scanpix.plot.multPlot`) are handled automatically: only the canvas is exported.

---

# B) Parameter space

_scaNpix_ uses three parameter spaces.

## 1. General (object) parameters

These are used when loading data and doing some basic pre-processing (e.g. position smoothing). You won't usually need to change most of them. They're stored in `obj.params` as a [map container](https://www.mathworks.com/help/matlab/map-containers.html). Use the name of a parameter as the key to get its value (e.g. `obj.params('posFs')` gives you the position sample rate), and `obj.params.keys` lists all parameters in the container.

The default values come from `scanpix.helpers.defaultParamsContainer(type)` and you should leave that file as it is. You can save your own version to a file instead: `obj.saveParams('container', 'YourFile')` writes the current map container to disk. Store your parameter file in `PathOnYourDisk\+scanpix\files\YourFile.mat`.

### Full list DACQ
* _scalePos2CamWin_: true/false. If true, rate maps are binned to the camera window (`xmin`/`xmax`/`ymin`/`ymax` from the `.set` file). If false (default) and `envSize` + `trialType` are known for a trial (e.g. from the optional meta xml file), positions are fitted to the environment instead.
* _ScalePos2PPM_: scale position data to this pix/m (_default=400_). This is particularly useful for keeping rate map sizes in proportion when you recorded in different environments with different sizes and/or pix/m settings for the camera.
* _posMaxSpeed_: speeds > posMaxSpeed are treated as tracking errors and ignored (set to _NaN_); in m/s (_default=4_)
* _posSmooth_: smooth position data over this many seconds (_default=0.4_)
* _posHead_: position of the head relative to the headstage LEDs (_default=0.5_)
* _cutFileType_: type of cut file, i.e. `'cut'` (Tint) or `'clu'` (KlustaKwik); _default='cut'_
* _cutTag1_: cut file tag that follows the base filename but precedes `_tetrodeN` in the filename (_default=''_). This is historic (and idiosyncratic) to the data of the original pup replay study, so you can probably ignore it.
* _cutTag2_: cut file tag that follows `_tetrodeN` in the filename (_default=''_)
* _APFs_: spike data sampling rate in Hz (_default=48000_)
* _lfpFs_: sampling rate of the `.eeg` files in Hz (_default=250_)
* _lfpHighFs_: sampling rate of the high sampling rate `.egf` files in Hz (_default=4800_)
* _loadHighFsLFP_: true/false; also load the high sampling rate EEG files (_default=true_)
* _loadAllWFs_: true/false. If true (default), every single waveform is kept. If false, only the mean waveform per cell is stored (saves a lot of memory).
* _defaultDir_: default directory where UI dialogues look for things, e.g. data (_default='Path/To/The/+scanpix/Code/On/Your/Disk'_)
* _myRateMapParams_: `'FileNameOfYourRateMapParams.mat'` (_default=''_)

### Full list Neuropixels data
* _scalePos2CamWin_: see DACQ (_default=false_).
* _scalePos2Env_: true/false. If true (and _scalePos2CamWin_ is false), positions are fitted to the physical size of the environment using `envSize` from the metadata xml (`scanpix.maps.scalePosition`), so rate maps have a fixed size per environment. If false (default), positions are only scaled to _ScalePos2PPM_ and rate maps span the visited area.
* _ScalePos2PPM_: scale position data to this pix/m (_default=400_). This is particularly useful for keeping rate map sizes in proportion when you recorded in different environments with different sizes and/or pix/m settings for the tracking.
* _posMaxSpeed_: speeds > posMaxSpeed are treated as tracking errors and ignored (set to _NaN_); in m/s (_default=4_)
* _posSmooth_: smooth position data over this many seconds (_default=0.4_)
* _maxPosInterpolate_: maximum duration of a gap of missing positions that is interpolated over; longer gaps are left as _NaN_; in s (_default=2.5_)
* _InterpPos2PosFs_: true/false; interpolate position data to the exact sampling rate (which is slightly different from exactly 50Hz). This substantially speeds up making rate maps (_default=true_).
* _posHead_: position of the head relative to the headstage LEDs (_default=0.5_)
* _posFs_: nominal position sampling rate in Hz (_default=50_)
* _loadFromPhy_: logical flag that sets which sorting results to use. If _true_ (default), the code tries Phy; otherwise it uses the raw Kilosort results.
* _APFs_: sampling rate of the Neuropixels AP band (_default=30000Hz_)
* _lfpFs_: sampling rate of the Neuropixels LFP band (_default=2500Hz_)
* _defaultDir_: default directory where UI dialogues look for things, e.g. data (_default='Path/To/The/+scanpix/Code/On/Your/Disk'_)
* _myRateMapParams_: `'FileNameOfYourRateMapParams.mat'` (_default=''_)

### Full list behavioural data
* _scalePos2Env_, _ScalePos2PPM_, _posMaxSpeed_, _posSmooth_, _maxPosInterpolate_, _posHead_, _posFs_, _defaultDir_, _myRateMapParams_: same as for Neuropixels data (_scalePos2Env_ default=false)

> **Note on older saved objects/param files:** `maxPosInterpolate` used to be a distance in cm (default 15–30). It is now a **duration in s**. If you load an old parameter file, check this value, because e.g. 30 is now read as 30 s.

## 2. Parameters for maps

These are stored as a scalar Matlab structure in `obj.mapParams` (a hidden property), with one substructure per map type (`obj.mapParams.rate`, `.dir`, `.lin`, `.linpos`, `.objVect`, `.speed`, `.sac`, `.gridProps`). The default values come from `scanpix.maps.defaultParamsRateMaps`. Again, don't edit anything in there.

You can change these parameters on the fly when making different kinds of maps (e.g. `obj.mapParams.rate.binSizeSpat = 2;` followed by `obj.addMaps('rate')`).

If you want your own values as the default, change them in the object and save them to disk with `obj.saveParams('maps', 'YourFile')` (they go to `PathOnYourDisk/+scanpix/files/YourFile.mat`). Then set `obj.params('myRateMapParams') = 'YourFile'`.

### Full list
* General params:
  * _speedFilterLimits_: limits for speed filtering in cm/s (_default=[2.5 400]_); copied into each map type below
  * _showWaitBar_: show a waitbar (_default=false_); copied into each map type below
* 2D rate maps (`.rate`):
  * _speedFilterFlagRMaps_: logical flag; speed filter the position data (_default=true_)
  * _speedFilterLimitLow_: lower speed limit in cm/s (_default=2.5_)
  * _speedFilterLimitHigh_: upper speed limit in cm/s (_default=400_)
  * _binSizeSpat_: bin size for spatial rate maps in cm (_default=2.5_)
  * _smooth_: type of smoothing; `'adaptive'` (default) or `'boxcar'`
  * _kernel_: size of the boxcar smoothing kernel in bins (_default=5_)
  * _alpha_: alpha parameter for adaptive smoothing (_default=200_; probably shouldn't be changed)
  * _trimNaNs_: trim rows or columns of the map that are all _NaN_ (_default=false_)
* Directional maps (`.dir`):
  * _speedFilterFlagDMaps_: logical flag; speed filter the position data (_default=true_)
  * _speedFilterLimitLow_ / _speedFilterLimitHigh_: speed limits in cm/s (_default=2.5 / 400_)
  * _binSizeDir_: bin size for directional maps in degrees (_default=6_)
  * _dirSmoothKern_: size of the smoothing kernel for directional maps in bins (_default=5_)
* Linear rate maps (`.lin`):
  * _binSizeLinMaps_: bin size for linear rate maps in cm (_default=2.5_)
  * _smoothFlagLinMaps_: logical flag; smooth the maps (_default=true_)
  * _smoothKernelSD_: SD of the Gaussian smoothing kernel in bins (_default=2_). The kernel is 5\*SD long.
  * _speedFilterFlagLMaps_: logical flag; speed filter the position data (_default=true_)
  * _speedFilterLimitLow_ / _speedFilterLimitHigh_: speed limits in cm/s (_default=2.5 / 400_)
  * _posIsCircular_: logical flag; treat position data as circular, e.g. on a square track (_default=false_)
  * _remTrackEnds_: set this many bins to NaN at each end of the linear track (_default=0_). Don't use for square track data.
* Linearisation of position data (`.linpos`):
  * _minDwellForEdge_: minimum dwell of the animal in a bin at the edge of the environment, in s (_default=1_)
  * _durThrCohRun_: minimum duration of a run in one direction, in s (_default=2_). Set to 0 if you don't want to remove position data from runs shorter than the threshold.
  * _filtSigmaForRunDir_: SD of the Gaussian filter that pre-filters the data before finding CW and CCW runs, in s (_default=3_). The kernel is 2\*SD long.
  * _durThrJump_: threshold for short changes of running direction, which are ignored if the gradient < gradThrForJumpSmooth; in s (_default=2_)
  * _gradThrForJumpSmooth_: gradient of the smoothed linear positions in the jump window, in cm/s (_default=2.5_)
  * _runDimension_: estimated from the data if left empty (_default=[]_). Only used for linear track data (a somewhat historic parameter, since it can usually be estimated from the data).
  * _dirTolerance_: tolerance for heading direction relative to the arm direction when calculating run direction on the track (_default=70°_)
* Object vector maps (`.objVect`):
  * _binSz_dist_: distance bin size in cm (_default=2.5_)
  * _binSz_dir_: direction bin size in degrees (_default=5_)
  * _smKernelSz_OV_: size of the smoothing kernel (_default=5_)
  * _smSigma_OV_: sigma of the smoothing kernel (_default=2_)
* Speed maps (`.speed`):
  * _minBinProp_: valid speed bins need to hold more than this proportion of all samples (_default=0.005_)
  * _binSizeSpeed_: speed bin size in cm/s (_default=2_)
  * _maxSpeed_: maximum speed in cm/s (_default=40_)
  * _confInt_: confidence interval in % (_default=95_)
* Spatial autocorrelograms (`.sac`):
  * _method_: method used to compute the autocorrelation; `'moser'` or `'barry'` (_default='moser'_)
  * _removeMinOverlap_: remove bins in the sAC where the rate map overlap was less than 20 bins (_default=true_)
  * _smooth_: smooth the sAC; not recommended (_default=false_)
  * _hSize_: size of the smoothing kernel (_default=5_)
  * _sigma_: sigma of the smoothing kernel (_default=1.5_)
* Grid properties (`.gridProps`):
  * _binAC_: bin the sAC (or not) before peak detection (_default=true_)
  * _nBinSteps_: number of bin steps (between -1 and 1) for binning the sAC (_default=21_)
  * _thresh_: set bins < _thresh_ to NaN (_default=0_)
  * _minPeakSz_: minimum number of pixels in a single peak (_default=8_)
  * _plotEllipse_: plot the ellipse fitted to the inner ring of the sAC (_default=false_)
  * _verbose_: switch on verbose mode (_default=false_)

## 3. Parameters for loading Neuropixels LFP

These are stored as a struct in `obj.lfpParams` (a hidden property). The defaults come from `scanpix.helpers.defaultParamsLFP`. Neuropixels LFP is large, so by default only a subset of channels is loaded:

* _chanSpacing_: spacing between loaded channels along the probe, in µm (_default=100_)
* _depthRange_: `'cells'` (default; span from the most dorsal to the most ventral cell in `obj.cell_ID(:,2)`, falling back to `'probe'` if no spikes are loaded yet), `'probe'` (whole probe) or `[minDepth maxDepth]` in µm from the probe tip
* _chans_: explicit list of probe channels (1-based) to load. Overrides _chanSpacing_ and _depthRange_ if not empty (_default=[]_).
* _downsampleFactor_: integer factor to downsample the LFP by (_default=1_, i.e. native 2.5kHz)

```matlab
obj.lfpParams.chanSpacing = 40;
obj.load({'lfp'});
```

---

# C) Class properties ( property (format) )

## General
* _params_ (containers.Map): map container with the general params
* _chanMap_ (struct): Kilosort channel map (Neuropixels only)
* _dataPath_ (string): `FullPathToDataOnDisk` of the raw data for each trial
* _dataPathSort_ (string): `FullPathToDataOnDisk` of the spike sorting output (usually the same as _dataPath_)
* _dataSetName_ (char): unique identifier for the dataset. For DACQ this is a nondescript name, because Axona data have no metadata file to get this information from.
* _trialNames_ (string array): list of filenames in the dataset
* _cell_ID_ (double): (nCells,4) array.
  * DACQ: cell ID (column 1), tetrode ID (column 2) and a simple numeric index (column 3)
  * Neuropixels: cluster ID (column 1), cluster depth (column 2) and the channel closest to the centre of mass of the cluster (column 3)
* _cell_Label_ (string): (nCells,1) string array with the label for each cluster from Kilosort (`'good'` or `'mua'`) or Phy (`'good'`, `'mua'` or `'noise'`) (Neuropixels only)
* _histo_reconstruct_ (struct): data for the histological reconstruction (Neuropixels only). One field per brain region, with the depth range, channels, cells and coverage.
* _trialMetaData_ (struct): non-scalar structure with trial-specific metadata (from e.g. `.set`, `.meta` or the xml files)
  * DACQ:
    * _tracked_spots_: number of LEDs
    * _xmin_, _xmax_, _ymin_, _ymax_: camera window
    * _sw_version_: software version
    * _trial_time_: start of the recording as time of day (as set on the machine)
    * _ADC_fullscale_mv_: scale max for channels at gain=1 in mV (for USB systems this should be 1.5V)
    * _lightBearing_: direction of the LEDs in degrees (up to 4 lights)
    * _colactive_: logical index of the active LEDs (probably only relevant for multi-colour LED tracking in DACQ)
    * _gains_: nTetrodes x 4 array of channel gains (up to 32 tetrodes (128 channels) possible)
    * _fullscale_: nTetrodes x 4 array of the scale max in µV
    * _lfp_channel_: nEEGs x 1 array of the channels the EEGs were recorded from
    * _lfp_recordingChannel_: nEEGs x 1 array of the channels set to EEG in DACQ (same as above if the EEG was recorded in mode SIGNAL, different if it was REF)
    * _lfp_slot_: nEEGs x 1 array of the EEG number in DACQ (so .eeg, .eeg2, … , .eegN)
    * _lfp_scalemax_: nEEGs x 1 array of the scale max for the EEG channels
    * _lfp_filter_: nEEGs x 1 array of the EEG filter type (0=DIRECT, 1=DIRECT+NOTCH, 2=HIGHPASS, 3=LOWPASS, 4=LOWPASS+NOTCH)
    * _lfp_filtresp_, _lfp_filtkind_, _lfp_filtfreq1_, _lfp_filtfreq2_, _lfp_filtripple_: settings of the user-defined filter (response type, kind (most likely Butterworth), lower/upper bounds, ripple)
    * _posFs_: position sampling rate
    * _ppm_: pixel/m. This holds the final ppm of the position data, so it differs from the original when you scale the data to a common ppm value.
    * _ppm_org_: pixel/m of the raw data
    * _PosIsScaled_: position is scaled to the standard ppm; yes/no
    * _PosIsFitToEnv_: position was fitted to the environment; `{flag, info}`
    * _duration_: duration of the trial in s
    * any extra fields from the optional `meta*_<trialName>.xml` (e.g. _envSize_, _trialType_, _envBorderCoords_)
    * _trackType_: `'sqtrack'` or `'linear'` (_default=''_)
    * _trackLength_: track length in cm (_default=[]_, since it differs for each type of track). For a square track, give the length of 1 arm only. This is crucial for matching rate map sizes across datasets.
  * Neuropixels data:
    * _log_: log for position data loading. It gives an overview of how good/not so good the tracking data are (missed sync pulses, missing frames, corrupt frame counts, interpolation, ...).
    * _animal_: animal number
    * _date_: date of data collection
    * _age_: age of the animal
    * _filename_: filename for the trial
    * _trialType_: identifier for the trial type (e.g. `'fam'` for familiar environment)
    * _duration_: duration of the trial in s
    * _envSize_: size of the recording environment in cm
    * _envBorderCoords_: coordinates of the recording environment borders in pix
    * _nLEDs_: number of LEDs used for tracking
    * _LEDfront_: which colour LED was at the front (if you used 2)
    * _LEDorientation_: orientation of the LEDs relative to the animal's head (if you used 2)
    * _objectPos_: xy coordinates of objects (if there were any in a trial)
    * _posFs_: actual position sampling rate
    * _ppm_: pixel/m. This holds the final ppm of the position data, so it differs from the original when you scale the data to a common ppm value.
    * _ppm_org_: pixel/m of the raw data
    * _nChanAP_ / _nChanSort_ / _nChanTot_: number of channels in the AP stream, used for sorting (usually 383), and in total (usually 385 for SpikeGLX or 384 for OE)
    * _missedSyncPulses_: list of missing pulses (ideally empty)
    * _offSet_: time of the first sync pulse in the raw data
    * _BonsaiCorruptFlag_: Bonsai data corrupt; yes/no (this can happen in various ways)
    * _PosIsScaled_: position is scaled to the standard ppm; yes/no
    * _PosIsFitToEnv_: position was fitted to the environment
    * _lfpFs_: sample rate of the stored LFP (after downsampling)
    * _lfpUVPerBit_: [1 x nChannels] conversion factor from raw LFP to µV
  * Behavioural data: similar fields to Neuropixels (from the xml + Bonsai file)
* _posData_ (struct): scalar structure with the position data. Fields:
  * _XYraw_: cell arrays of raw LED position data in pixels (xy coordinates)
  * _XY_: cell arrays of the processed animal position in pixels (xy coordinates)
  * _direction_: cell arrays of the animal's head direction in degrees. All data types use the same convention: y-axis pointing down, i.e. angles run clockwise on screen, so directional maps line up with rate maps. **Note:** DACQ data loaded with older versions used the Tint (anticlockwise) convention.
  * _speed_: cell arrays of the animal's running speed in cm/s
  * _linXY_: cell arrays of linearised position (only made when you make linear rate maps)
  * _sampleT_: sample times of the position frames grabbed from the Bonsai data. Not really used for anything (Neuropixels only).
* _spikeData_ (struct): scalar structure with the spike data. Fields:
  * _spk_Times_: cell arrays of spike times (in s)
  * _spk_waveforms_: cell arrays of waveforms (in µV). The format for each cell is nSpikes x nSamples x nChannels (i.e. nSpikes x 50 x 4 for DACQ). With `loadAllWFs=false` (DACQ), only the mean waveform is kept.
  * _sampleT_: timestamps of the position frames in Neuropixels time. Only relevant for Neuropixels data.
* _lfpData_ (struct): scalar structure with the LFP data. Fields:
  * _lfp_: cell arrays of LFP.
    * DACQ: low sample rate (250Hz) EEG in µV
    * Neuropixels: **raw `int16`** [nChannels x nSamples] to save memory. Use `scanpix.lfp.lfp2uV` (or multiply by `trialMetaData.lfpUVPerBit`) to convert to µV. Sample 1 is the sample closest to the start of the trial.
  * _lfpHighSamp_: cell arrays of high sample rate (4800Hz) EEG in µV (DACQ only)
  * _lfpTet_: cell array of the tetrode IDs the EEGs were recorded from (DACQ only)
  * _lfpChans_: [probe channel (1-based), depth from tip in µm] for each row of _lfp_, sorted by depth (Neuropixels only)
* _bhaveData_ (struct): behavioural data (behavioural objects only)
* _maps_ (struct): scalar structure with the maps
  * _rate_: cell arrays of standard 2D rate maps
  * _spike_: cell arrays of 2D spike maps
  * _pos_: cell arrays of 2D position maps
  * _dir_: cell arrays of directional maps
  * _sACs_: cell arrays of spatial autocorrelograms
  * _OV_: cell arrays of object vector maps
  * _speed_: cell arrays of speed maps
  * _lin_: cell(3,1) arrays of linearised rate maps: the full rate map ({1}), the rate map for CW runs ({2}) and for CCW runs ({3})
  * _linPos_: cell arrays of linearised position maps, in the same format as _lin_

## Hidden properties
* _fileType_ (char): identifier for the data file type (`.set` for DACQ, `.ap.bin` for Neuropixels and `.csv` for behavioural data)
* _type_ (char): data type of the object (`'dacq'`, `'npix'`, `'bhave'`, `'nexus'`)
* _fields2spare_ (cell array): fields that aren't changed when e.g. deleting or reordering data in the object. Typically these are fields with only 1 value per dataset (_default={'params','dataSetName','cell_ID','cell_Label','histo_reconstruct'}_). If you add new properties that should be spared, add them here!
* _mapParams_ (struct): map params (see `scanpix.maps.defaultParamsRateMaps`)
* _lfpParams_ (struct): params for loading Neuropixels LFP (see `scanpix.helpers.defaultParamsLFP`)
* _loadFlag_ (logical): flag that shows whether any data has been loaded into the object (_default=false_)
* _isConcat_ (logical): flag that shows whether the data were loaded from a concatenated file (_default=false_)

   
 



