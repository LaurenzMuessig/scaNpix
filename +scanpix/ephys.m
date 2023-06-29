classdef ephys < handle
    % ephys data class. Let's you load in different types of ephys data into a class
    % object for analysis or data inspection
    %
    % Class constructor: scanpix.ephys
    % Construct class obj
    %
    % Usage:
    %       obj = scanpix.ephys;
    %       obj = scanpix.ephys(type);
    %       obj = scanpix.ephys(type,prmsMode);
    %       obj = scanpix.ephys(type,prmsMode,setDirFlag);
    %
    % Inputs:
    %       type:       data type ('dacq', 'npix', 'nexus')
    %
    %       prmsMode:   'default' - uses default parameter (default)
    %                   'ui'      - opens UI dialogue
    %                   'file'    - load params from file
    %
    %       setDirFlag: true (default) or false
    %                   - if false will skip UI dialogue for data selection
    %                   (e.g. when you use constructor programatically)
    %
    % LM 2020
    
    
    %% PROPERTIES
    properties % params %
        % Map object
        params                containers.Map
        chanMap               struct
    end
    
    properties % meta data %
        dataPath(1,:)         string
        dataPathSort(1,:)     string
        dataSetName(1,:)      char
        trialNames(1,:)       string
        cell_ID(:,4)          double %{mustBePositive, mustBeNonNan, mustBeNonzero}
        cell_Label(:,1)       string
    end
    
    properties % trial data %
        trialMetaData(1,:)    struct
        posData               struct  = struct('XYraw',[],'XY',[],'direction',[],'speed',[],'linXY',[],'sampleT',[]);
        spikeData             struct  = struct('spk_Times',[],'spk_waveforms',[],'sampleT',[]); 
        lfpData               struct  = struct('lfp',[]);
    end
    
    properties % maps %
        maps                  struct  = struct('rate',[],'spike',[],'pos',[],'dir',[],'sACs',[],'OV',[],'speed',[],'lin',[],'linPos',[]);
    end
%     
    properties(Hidden)
        fileType              char    = '';
        type                  char {mustBeMember(type,{'npix','dacq','nexus','init'})} = 'init';
        %         uniqueCellIndex(:,1)  logical
        fields2spare          cell    = {'params','dataSetName','cell_ID','cell_Label'}; % spare this when deleting or rearranging data. make sure to add new properties that should be spared here!
        mapParams             struct  = scanpix.maps.defaultParamsRateMaps;
        loadFlag              logical = false;                                         % flag so we know something has been loaded into object
        badChans
    end
    
    %% METHODS
    % constructor %
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = ephys(type,prmsMode,setDirFlag)
            % dacq - creates class object
            %
            % Syntax:
            %       obj = npix;
            %       obj = npix(prmsMode);
            %       obj = npix(prmsMode, uiFlag);
            %
            % Inputs:
            %    prmsMode    - 'default' - uses default parameter (default)
            %                  'ui'      - opens UI dialogue
            %                  'file'    - load params from file
            %    setDirFlag  - true (default)/false - skip UI set file
            %                selection dialogue if false
            %
            % Outputs:
            %    obj      - ephys object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj = obj@handle;
            
            if nargin == 0
                % select data type for object
                str = {'dacq','npix','nexus'};
                [select, loadCheck] = listdlg('PromptString','Select what data tickles your fancy:','ListString',str,'ListSize',[160 100],'SelectionMode','Single');
                if ~loadCheck
                    warning('scaNpix::ephys:: Object creation aborted. And I thought you loved ephys class objects...');
                    return;
                else
                    obj.type = str{select};
                end
            else
                obj.type = type;
            end

            %
            switch obj.type
                case {'dacq','nexus'}
                    obj.fileType = '.set';
                case 'npix'
                    obj.fileType = '.ap.bin';
            end
            %
            if nargin <= 1
                obj.getParams('default');
                warning('scanNpix: No params selected. Using defaults.');
            else
                obj.getParams(prmsMode); 
            end
            if isempty(obj.params); return; end
            
            if nargin <= 2
                setDirFlag = true;
            end
            
            if setDirFlag
                obj.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
                if isempty(obj.trialNames); return; end
                obj.selectTrials(obj.trialNames);
            end
                       
            if ~isempty(obj.dataPath)
                obj.params('defaultDir') = [fileparts(obj.dataPath{1}) filesep];  % use parent directory rather than absolut as this makes more sense for multiple trial data
            end
            
        end
    end
    
    %%
    methods
        
        %%
        function changeParams(obj, changeMode, optionalPath2file)
            % changeParams - change params container in obj
            %
            % Syntax:
            %       obj.changeParams
            %       obj.changeParams(changeMode)
            %       obj.changeParams(changeMode, optionalPath2file)
            %
            % Inputs:
            %    changeMode - 'default' - uses default parameters;
            %                 'ui' - opens UI dialogue (default);
            %                 'file' - load params from file
            %                 (Note that this gets passed on as 'mode' to scanpix.helpers.getParams)
            %
            %
            %    optionalPath2file - string; optional name of .mat file with parameters (we assume it's in 'YourPathOnDisk\+scanpix\files\')
            %                        if ommited and changeMode='file' will open UI dialogue to fetch file
            %
            % Outputs:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin < 2
                changeMode = 'ui';
            end
            
            if nargin < 3
                optionalPath2file = [];
            end
            
            % change params
            obj.getParams(changeMode, optionalPath2file);
%             % in case we aborted loading we want to keep things as they are
%             % and don't change anything
%             if ~isempty(prmsMap)
%                 obj.params = prmsMap;
%             end
        end
        
        %%
        function saveParams(obj, type, fileName)
            % saveParams - saves current params map container or rate map params to disk as .mat file
            % User can save her/his own standard parameters
            %
            % Syntax:
            %       obj.saveParams(type)
            %       obj.saveParams(type, fileName)
            %
            % Inputs:
            %    type     - string; 'container' - save general params; 'maps' - save rate map params;
            %               Note files should be saved to 'PathOnDisk\+scanpix\files\' and that 'defaultParams' & 'defaultRatemapParams' are
            %               not allowed as filenames
            %
            %    fileName - string; filename for params.mat file
            %
            % Outputs:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            classFolder = what('scanpix');
            
            switch lower(type)
                case 'container'
                    prmsMap = obj.params; %#okagrow
                    outStr  = 'prmsMap';
                    suggest = 'MyParams';
                case 'maps'
                    prms    = obj.mapParams; %#okagrow
                    outStr  = 'prms';
                    suggest = 'MyRateMapParams';
                otherwise
                    ME = MException('scaNpix:ephys::saveParams:InvalidType', ['''' type ''' is not a valid type of parameter group you can save. Try ''container'' or ''maps'' instead']);
                    throw(ME);
            end
            
            if nargin < 3 || isempty(fileName)
                fileName     = {'defaultParams'};
                while any( strcmpi(fileName,{'defaultparams','defaultratemapparams'}) )
                    fileName = inputdlg('Please enter filename without extension (''defaultParams'' or ''defaultRateMapParams'' are not allowed as filenames)', 'Save Params', 1, {suggest} );
                    if isempty(fileName)
                        warning('scaNpix: saving params map container aborted');
                        return
                    end
                end
            else
                if any( strcmpi(fileName,{'defaultparams','defaultratemapparams'}) )
                    ME = MException('scaNpix:ephys::saveParams:InvalidFilename', '''defaultParams'' or ''defaultRateMapParams'' aren''t allowed as file names mate!');
                    throw(ME);
                end
            end
            
            save( fullfile(classFolder.path,'files',[ fileName{:} '.mat' ]), outStr);
        end
        
        %%
        function addMetaData(obj,name,values)
            % addMetaData - add meta data tags to object. These will be added to obj.trialMetaData.(name)
            %
            % Syntax:
            %       obj.addMetaData
            %       obj.addMetaData(name)
            %       obj.addMetaData(name, varargin)
            %
            % Inputs:
            %    name     - string, name of field to be added to obj.metaData. If ommited will open UI dialogue to fetch data
            %    values   - values for field 'name'. if single value is given for dataset containing multiple trials value will
            %               be expanded for all trials.
            % Outputs:
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin < 2
                % add meta data per UI input
                uiInput = inputdlg(['Metadata field name'; strcat(obj.trialNames,' - value')'],'Please indicate meta data field you want to add',1);
                if isempty(uiInput)
                    warning('scaNpix::ephys::addMetaData: Cancelled adding meta data to object. It would have been so nice to get some more of it...');
                    return;
                end
                name   = uiInput{1};
                values = uiInput(2:end);
                % convert numeric input to numbers
                for i = 1:length(values)
                    if any(ismember(values{i}, '0123456789+-.eEdD'))
                        values{i} = str2num(values{i}); %#okAgrow
                    end
                end
            end
            
            if nargin == 2
                % initialise empty field
                values = repmat({[]},1,length(obj.trialNames));
            end
            
            if ~iscell(values)
                values = {values};
            end
            
            if length(values) < length(obj.trialNames)
                values = repmat(values,1,length(obj.trialNames)); % expand
            end
            
            % add metadata
            for i = 1:length(obj.trialNames)
                obj.trialMetaData(i).(name) = values{i};
            end
        end
        
        %%
        function addData(obj,newTrialNames,loadMode,reorderFlag)
            % addData - add new data to object, i.e. load additional trials into dataset
            % (These will be added to end of sequence in object by default)
            %
            % Syntax:
            %       obj.addData
            %       obj.addData(newTrialNames)
            %       obj.addData(newTrialNames, loadMode)
            %       obj.addData(newTrialNames, loadMode, reorderFlag )
            %       obj.addData([],__)
            %
            % Inputs:
            %    newTrialNames - string(s) of set file to be added to object
            %                    if ommited or empty will open UI dialogue to fetch files
            %    loadMode      - cell array; {'all'} (default) or any combination of {'pos','spikes','lfp'}
            %    reorderFlag   - logical (default=true); if true run also obj.reorderData
            %
            %
            % Outputs:
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            warning('scaNpix::ephys::addData:All existing rate maps and metadata fields will be removed when adding data to object. You''ll need to re-generate these.'); % could add UI confirmation box?
            
            if ~obj.loadFlag
                warning('scaNpix::ephys::addData:You need to load data into the object first before you can use ''obj.addData''.');
                return;
            end
            
            if nargin < 2 || isempty(newTrialNames)
                [newTrialNames, newDataDir] = obj.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
                if isempty(newTrialNames)
                    return
                end
                [newTrialNames, ind] = obj.selectTrials(newTrialNames);
            else
                % this only works for 1 trial currently! 
                [path,newTrialNames] = fileparts(newTrialNames);
                newDataDir = {[path filesep]};
                ind = 1;
                if strcmp(obj.type,'npix')
                    [~,newTrialNames] = fileparts(newTrialNames);
                end
            end
            
            if nargin < 3
                loadMode = {'meta'};
                if ~isempty(obj.posData.XYraw{1})
                    loadMode = [loadMode(:); 'pos']; 
                end
                if ~isempty(obj.spikeData.spk_Times{1})
                    loadMode = [loadMode(:); 'spikes']; 
                end
                if ~isempty(obj.lfpData.lfp{1})
                    loadMode = [loadMode(:); 'lfp']; 
                end
            end
            
            if nargin < 4
                reorderFlag = true;
            end
            
            if any(strcmp(obj.trialNames,newTrialNames))
                warning('scaNpix::ephys::addData:The dataset you want to add is already in object, so it will just be reloaded.');
            else
                obj.trialNames = horzcat(obj.trialNames, newTrialNames);
                obj.dataPath   = horzcat(obj.dataPath, newDataDir(ind));
            end
            %
            obj.load(loadMode, newTrialNames);
            
            % reorder
            if reorderFlag
                obj.reorderData;
            end 
        end
        
        %%
        function deleteData(obj, type, varargin)
            % delete - delete data of individual trials from data objects.
            %
            % Syntax:
            %       obj.delete(type)
            %       obj.delete('trials',trialStr)
            %       obj.delete('cells',cellInd)
            %
            % Inputs:
            %    type     - 'trials' or 'cells'
            %    varargin:
            %             - trialStr - string/cell array of strings; name(s) of trial(s) to be deleted from object
            %               if ommited will open UI dialogue to select trial(s)
            %             - cellInd - logical or numeric index of cells to be deleted from data
            %
            % Outputs:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin == 1
                str = {'trials','cells'};
                [select, loadCheck] = listdlg('PromptString','Select what data you hate:','ListString',str,'ListSize',[160 100],'SelectionMode','Single');
                if ~loadCheck
                    warning('scaNpix::ephys::deleteData: Agreed. Best to keep your data where it is.');
                    return;
                else
                    type = str{select};
                end
            end
            
            % parse inputs
            switch type
                case 'trials'
                    if nargin < 3
                        [select, loadCheck] = listdlg('PromptString','Select what data to delete:','ListString',obj.trialNames,'ListSize',[160 100]);
                        if ~loadCheck
                            warning('scaNpix::ephys::deleteData: Deleting data aborted. More data is better anyway.');
                            return;
                        end
                        trialStr = obj.trialNames(select);
                    else
                        trialStr = varargin{1};
                    end
                case 'cells'
                    if nargin < 3
                        %             [select, loadCheck] = listdlg('PromptString','Select what data to delete:','ListString',obj.trialNames,'ListSize',[160 100]);
                        %             if ~loadCheck
                        %                 warning('scaNpix: Deleting data aborted. More data is better anyway.');
                        %                 return;
                        %             end
                        %             trialStr = obj.trialNames(select);
                    else
                        if islogical(varargin{1})
                            cellInd = ~varargin{1};
                        elseif isnumeric(varargin{1})
                            cellInd = true(size(obj.cell_ID,1),1);
                            cellInd([varargin{1}]) = false;
                        end
                    end
            end
            
            % delete
            switch type
                case 'trials'
                    deleteInd = ismember(obj.trialNames,trialStr);
                    fields = fieldnames(obj);
                    for i = 1:length(fields)
                        if ~any(strcmp(fields{i},obj.fields2spare)) && ~isempty(obj.(fields{i}))
                            if isstruct(obj.(fields{i})) && isscalar(obj.(fields{i}))
                                obj.(fields{i}) = structfun(@(x) x(~deleteInd),obj.(fields{i}),'uni',0);
                            else
                                obj.(fields{i}) = obj.(fields{i})(~deleteInd);
                            end
                        end
                    end
                    
                    % remove the last bits if object is now empty
                    if isempty(obj.trialNames)
                        for i = 1:length(obj.fields2spare)
                            if ~strcmp(obj.fields2spare{i},'params')
                                obj.(obj.fields2spare{i}) = [];
                            end
                        end
                        obj.loadFlag  = false;
                        obj.changeParams('default'); % reset defaults
                        obj.mapParams = scanpix.maps.defaultParamsRateMaps(class(obj)); % reset defaults
                        warning('scaNpix::ephys::deleteData: Back to square one...');
                    end
                    
                case 'cells'
                    
                    f = fieldnames(obj.spikeData);
                    for i = 1:length(f)
                        if ~isempty(obj.spikeData.(f{i}){1}) && ~strcmp(f{i},'sampleT')
                            obj.spikeData.(f{i}) = cellfun(@(x) x(cellInd),obj.spikeData.(f{i}),'uni',0);
                        end
                    end
                    % need to test whether this works for 2D rate maps as well
                    f = fieldnames(obj.maps);
                    for i = 1:length(f)
                        indEmpty = cellfun('isempty',obj.maps.(f{i}));
                        if ~all(indEmpty)
                            if all(cellfun(@(x) size(x,1),obj.maps.(f{i})(~indEmpty)) == length(cellInd)) % standard and lin pos maps should be omitted
                                obj.maps.(f{i})(~indEmpty) = cellfun(@(x) x(cellInd,:),obj.maps.(f{i})(~indEmpty),'uni',0);
                            end
                            % lin maps are a pain to deal with - there must be an easier solution!
                            if strcmp(f{i},'lin')
                                numInd = find(~indEmpty);
                                for j = 1:length(numInd)
                                    tmpMaps = cell(3,1);
                                    for k = 1:3
                                        tmpMaps(k,1)  = cellfun(@(x) x{k}(cellInd,:),obj.maps.(f{i})(numInd(j)),'uni',0);
                                    end
                                    obj.maps.(f{i})(numInd(j)) = {tmpMaps};
                                end
                            end
                        end
                    end
                    obj.cell_ID = obj.cell_ID(cellInd,:);
                    if strcmp(obj.type,'npix')
                        obj.cell_Label = obj.cell_Label(cellInd,:);
                    end
            end
        end
        
        %%
        function reorderData(obj, orderIndex)
            % reorderData - reorder trial sequence in data object e.g. after adding or removing data
            %
            % Syntax:
            %       obj.reorderData
            %       obj.reorderData(orderIndex)
            %
            % Inputs:
            %    orderIndex - numeric index; will order data in object accordingly
            %                 (if ommited will open UI dialogue to indicate index)
            %
            %
            % Outputs:
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin < 2
                uiInput = inputdlg({obj.trialNames{:}},'Please indicate new index for each trial');
                if isempty(uiInput)
                    warning('scaNpix::ephys::reorderData: Reordering of trials aborted. Just ask yourself why you started it then...');
                    return;
                end
                orderIndex = str2num(cell2mat(uiInput))';
                
                if any(orderIndex > length(obj.trialNames))
                    ME = MException('scaNpix:ephys::reorderData:orderIndex','OrderIndex for trials > number of trials. That can''t work, can it?');
                    throw(ME);
                end
            end
            
            % re-order
            fields = fieldnames(obj);
            for i = 1:length(fields)
                if ~any(strcmp(fields{i},obj.fields2spare)) && ~isempty(obj.(fields{i}))
                    if isstruct(obj.(fields{i})) && isscalar(obj.(fields{i}))
                        obj.(fields{i}) = structfun(@(x) x(orderIndex),obj.(fields{i}),'uni',0);
                    else
                        obj.(fields{i}) = obj.(fields{i})(orderIndex);
                    end
                end
            end
            
        end
        
    end
    
    %% data loading %%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function load(obj, loadMode, varargin)
            % load - load a dataset(s) into class object.
            %
            % Syntax:
            %       obj.load
            %       obj.load(loadMode)
            %       obj.load(loadMode, dataset1, dataset2,...,datasetN)
            %
            % Inputs:
            %    loadMode - cell array; {'all'} or any combination of {'pos','spikes','eeg'}
            %    varargin - list of set file names (ommit extensions)
            %
            % Outputs:
            %
            % see also: obj.loadMetaData; obj.loadPosData; obj.loadSpikeData; obj.loadLFPData
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            reloadFlag = false;
            check0spikeCells = false;
            
            if strcmp(obj.type,'init')
                warning('scaNpix::ephys::load: You haven''t set the data type for the object yet (and probably neither the param space. Are you serious?');
            end
            
            if isempty(obj.dataPath)
                obj.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
                if isempty(obj.trialNames); return; end
                obj.selectTrials(obj.trialNames);
            end
            
            if nargin < 2
                str = {'all','pos','spikes','lfp'};
                [select, loadCheck] = listdlg('PromptString','Select what data to load:','ListString',str,'ListSize',[160 100]);
                if ~loadCheck
                    warning('scaNpix::ephys::load: No data selected. Nothing is loaded. Boring...');
                    return;
                else
                    loadMode = str(select);
                end
            end
            
            if ~iscell(loadMode)
                loadMode  = {loadMode};
            end
            
            if any(strcmp(loadMode,'reload'))
                reloadFlag = true;
                loadMode = loadMode(~strcmp(loadMode,'reload'));
            end
            
            if nargin < 3
                loadStr = obj.trialNames;
            else
                loadStr = obj.trialNames( ismember(obj.trialNames,varargin{1}) );
            end
            
            trialInd = find(ismember(obj.trialNames,loadStr));
            
            for i = 1:length(loadStr)
                
                % load some meta data
                if ~reloadFlag; obj.loadMetaData(trialInd(i)); end
                               
                for j = 1:length(loadMode)
                    
                    switch lower(loadMode{j})
                        case 'all'
                            obj.loadPositionData(trialInd(i));
                            obj.loadSpikeData(trialInd(i));
                            %                             obj.loadLFPData(trialInd(i));
                            break;
                        case 'pos'
                            obj.loadPositionData(trialInd(i));
                        case 'spikes'
                            if strcmp(obj.type,'dacq')
                                cellTet{i,1}     = obj.loadSpikeData(trialInd(i));
                                check0spikeCells = true; % flag
                            else
                                obj.loadSpikeData(trialInd(i),reloadFlag);
                            end
                        case 'lfp'
                            obj.loadLFPData(trialInd(i)); 
                        case 'meta'
                            % already loaded!
                        otherwise
                            ME = MException('scaNpix:ephys::load:invalidInput', ['' loadMode{j} ''' is not a valid data type selector. Try ''all'', ''pos'', ''spikes'' and/or ''eeg'' instead. ']);
                            throw(ME);
                    end
                end
            end
            
            % make sure to add empty cells to all trials if a unit fires 
            % only on some trials
            if check0spikeCells
                cellTet  = vertcat( cellTet,{ obj.cell_ID(:,1:2) } );
                allCellN = sortrows(unique(vertcat(cellTet{:}),'rows'),2);
                
                for i = 1:sum(~cellfun('isempty',cellTet))
                    % make sure that missing cells get added to spiketime cell arrays
                    missCellInd                        = ~ismember(allCellN,cellTet{i},'rows');
                    if ~any(missCellInd), continue, end
                    % add empty cells into trial data if necessary
                    temp                               = cell(size(allCellN,1), 1);
                    temp(~missCellInd)                 = obj.spikeData.spk_waveforms{i};
                    obj.spikeData.spk_waveforms{i}     = temp; %final output
                    temp(~missCellInd)                 = obj.spikeData.spk_Times{i};
                    obj.spikeData.spk_Times{i}         = temp; %final output
                end
                obj.cell_ID = [ allCellN ( 1:size(allCellN,1) )' nan(size(allCellN,1),1)]; % also add a numerical index for cells (which is sometimes useful to have)
            end

            obj.loadFlag = true; % flag loading
            
            % pre-allocate some data (e.g. part that isn't loaded), such
            % that all properties have same size
            obj.preallocEmpty(true,{'posData','spikeData','lfpData'});
        end
    end
    
    methods(Hidden)
        %%
        function loadMetaData(obj, trialIterator)
            % loadMeta - load meta data
            % Method for ephys class objects (hidden)
            %
            % Syntax:  obj.loadMetaData(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch obj.type
                case 'dacq'
                    scanpix.dacqUtils.loadSet(obj,trialIterator);
                case 'npix'
                    scanpix.npixUtils.loadMeta(obj,trialIterator);
                case 'nexus'
            end
        end
        
        %%
        function loadPositionData(obj, trialIterator)
            % loadPos - load position data
            % Method for ephys class objects (hidden)
            %
            % Syntax:  obj.loadPos(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch obj.type
                case 'dacq'
                    scanpix.dacqUtils.loadPos(obj,trialIterator);
                case 'npix'
                    scanpix.npixUtils.loadPosNPix(obj,trialIterator);
                case 'nexus'
            end
        end
        
        %%
        function cellIDs = loadSpikeData(obj, trialIterator, reloadFlag)
            % loadSpikes - load spike data
            % Method for ephys class objects (hidden)
            %
            % Syntax:  obj.loadSpikes(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% To-Do
            switch obj.type
                case 'dacq'
                    cellIDs = scanpix.dacqUtils.loadSpikes(obj,trialIterator);
                case 'npix'
%                     scanpix.npixUtils.loadSpikesNPix(obj,trialIterator,reloadFlag, false);
                    scanpix.npixUtils.loadSpikesNPix(obj,trialIterator,reloadFlag);
                case 'nexus'
            end
            
        end
        
        %%
        function loadLFPData(obj, trialIterator)
            % loadSpikes - load spike data
            % Method for ephys class objects (hidden)
            %
            % Syntax:  obj.loadSpikes(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% To-Do
            switch obj.type
                case 'dacq'
                    scanpix.dacqUtils.loadLFPs(obj,trialIterator);
                case 'npix'
                    %%%%
                case 'nexus'
            end
            
        end
        
        %%
        function prmsMap = getParams(obj,mode,path2file)
            % getParams - grab params container for ephys class object
            %
            % Syntax:
            %       prmsMap = obj.getParams
            %       prmsMap = obj.getParams(mode)
            %       prmsMap = obj.getParams(mode, path2file)
            %
            % Inputs:
            %    mode      - 'default' - uses default parameters (default)
            %              - 'ui'      - opens UI dialogue to grab parameters
            %              - 'file'    - load params from file
            %    path2file - optional name of .mat file with parameters (we assume it's in 'PathOnDisk\+scanpix\files\')
            %                if ommited and mode='file' will open UI dialogue to fetch file
            %
            % Outputs:
            %    prmsMap   - params map container
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin == 1 || strcmpi(mode,'default')
                
                prmsMap = scanpix.helpers.defaultParamsContainer(obj.type);
                
            elseif strcmpi(mode,'ui')
                % gather input from user
                if isempty(obj.params)
                    prmsMap = scanpix.helpers.defaultParamsContainer(obj.type); % we'll use the default params as base and then create a dialogue from these
                else
                    prmsMap = obj.params; % we'll use the current ones as base in case user wants to update existing ones
                end
                prompts     = prmsMap.keys;
                defaultVals = prmsMap.values;
                output      = scanpix.helpers.makeCustomUIDialogue(prompts, defaultVals);
                % exit gracefully
                if isempty(output)
%                     prmsMap = containers.Map;
                    return
                end
                % create the container
                prmsMap  = containers.Map(output(:,1),output(:,2));
                
            elseif strcmpi(mode,'file')
                % load from file
                if nargin < 3 || isempty(path2file)
                    classFolder = what('scanpix');
                    [fNames, dataDir] = uigetfile(fullfile(classFolder.path,'files','*.mat'),'Select params containers Map to load.');
                    if isnumeric(fNames)
                        warning('scaNpix::ephys::getParams: Loading of params container aborted!');
%                         prmsMap = containers.Map;
                        return;
                    end
                    path2file         = fullfile(dataDir,fNames);
                end
                temp                  = load(path2file);
                f                     = fieldnames(temp);
                prmsMap               = temp.(f{1});
                
            else
                ME = MException('scaNpix:ephys::getParams:invalidMode', [mode ' is not an accepted input to construct parameter container. Try ''default'',''ui'' or ''file'' instead.']);
                throw(ME);
            end
            
            if ~strcmp(prmsMap('type'),obj.type)
                ME = MException('scaNpix:ephys::getParams:invalidPRMSMap', [prmsMap('type') ' and ' obj.type ' don''t match. I knew this would happen but you didn''t want to listen.']);
                throw(ME);
            end

            obj.params = prmsMap;
        end
        
        %%
        function [trialNames, dataDir] = fetchFileNamesAndPath(obj,fileExt, defaultDir)
            % fetchFileNamesAndPath - grab filenames for trials and path on disk for data
            % The idea is that we choose a top level directory and then we'll find all
            % data files with '.fileExt' in the subdirectories. We will use a UI dialogue
            % to fetch that.
            %
            % Syntax:
            %    obj.fetchFileNamesAndPath(fileExt)
            %    obj.fetchFileNamesAndPath(fileExt, defaultDir )
            %
            % Inputs:
            %    defaultDir - path to directory where UI dialogue will start
            %
            % Outputs:
            %    trialNames - cell array of data file names (w/o extension)
            %    dataDir    - full path(s) to data on disk
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin == 1
                classFolder = what('scanpix');
                defaultDir = [classFolder.path filesep];
            end
            
            dirList = scanpix.fxchange.uigetfile_n_dir(defaultDir,'Select Top Level Folder(s) With Session Data');
            if isempty(dirList)
                warning('scaNpix::ephys::fetchFileNamesAndPath: Loading Aborted. Und tschuess!');
                [trialNames, dataDir] = deal(NaN);
                return
            end
            
            % now gather all datasets - search all subfolders for data files
            trialStruct = struct(dir('1'));
            for i = 1:length(dirList)
                trialStruct = vertcat(trialStruct,dir(fullfile(dirList{i},'**',['*' fileExt])));
            end
            
            if isempty(trialStruct)
                warning('scaNpix::ephys::fetchFileNamesAndPath: No ''%s'' file(s) found in selcted folder(s). Go and find your data yourself.',fileExt);
                return
            end
            
            [~,ind]      = sort([trialStruct.datenum]);
            tempNames    = {trialStruct(ind).name}';
            C = regexp(tempNames,fileExt,'split');
            
            if nargout > 1
                dataDir = strcat({trialStruct(ind).folder}', filesep);
                % remove extentions
                trialNames = cellfun(@(x) x{1}, C,'uni',0);
            else
                obj.dataPath = strcat({trialStruct(ind).folder}', filesep);
                % remove extentions
                obj.trialNames = cellfun(@(x) x{1}, C,'uni',0);
            end
            
            
        end
        %%
        function [trialNameStrOut, select] = selectTrials(obj,trialNameStrIn)
            % selectTrials - select from a list which trials to choose for loading data into ephys class objects
            %
            % Syntax:
            %       obj.selectTrials(trialNameStrIn)
            %
            % Inputs:
            %    trialNameStrIn  - char/cell array of filename strings that we want to
            %                      choose from
            %
            % Outputs:
            %    trialNameStrOut - user choice of 'trialNameStrIn'
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % skip if only one trial in object
            if ischar(trialNameStrIn)
                trialNameStrOut = trialNameStrIn;
                select = true;
                return;
            end
            
            % UI selection
            [select, loadCheck] = listdlg('PromptString','Select which Trial(s) to Include:','ListString',trialNameStrIn,'ListSize',[200 400],'CancelString','Keep All');
            if ~loadCheck
                return;
            else
                trialNameStrOut = trialNameStrIn(select);
                
                if nargout == 1
                    ind            = ismember(obj.trialNames,trialNameStrOut);
                    obj.trialNames = obj.trialNames(ind);
                    obj.dataPath   = obj.dataPath(ind);
                end
            end
            
        end
        %%
        function preallocEmpty(obj,noLoadProps,loadProps)
            % preallocEmpty - preallocate empty properties in ephys class object
            % (e.g. from things that weren't loaded) that have a trial based structure
            % (i.e. all raw data like positions, spiketimes etc.)
            %
            % Syntax:
            %       obj.preallocEmpty(noLoadProps)
            %       obj.preallocEmpty(noLoadProps,loadProps)
            %
            % Inputs:
            %    noLoadProps - true/false - flag to do properties that won't be loaded
            %                  through normal routines
            %    loadProps   - cell array; any combination of {'pos','spikes','lfp'}
            %
            % Outputs:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin < 3
                loadProps   = {};
            end
            
            % these will never get loaded automatically by normal loading routine
            if noLoadProps
                obj.posData.linXY = cell(1,length(obj.trialNames));
                obj.maps          = structfun(@(x) cell(1,length(obj.trialNames)),obj.maps,'uni',0);
            end
            
            % this could actually have a hard coded list instead of supplying as input?
            if ~isempty(loadProps)
                for i = 1:length(loadProps)
                    f = fieldnames(obj.(loadProps{i}));
                    for j = 1:length(f)
                        if isempty(obj.(loadProps{i}).(f{j})) || isempty(obj.(loadProps{i}).(f{j}){1})
                            obj.(loadProps{i}).(f{j}) = cell(1,length(obj.trialNames));
                        end
                    end
                end
            end
        end
    end
    
    methods(Static)
        
%         function dependentInd = checkPropDepend
%             metaClass = ?scanpix.ephys;
%             dependentInd = [metaClass.PropertyList.Dependent];
%             dependentInd = dependentInd(~[metaClass.PropertyList.Hidden]); 
%         end
    end
    

    
    methods
        
        %%
        function addMaps(obj, mapType, trialInd, varargin )
            % addMaps - create different types of firing rate maps for data in object
            %
            % We can add spatial rate maps ('rate'), directional rate maps ('dir') or
            % linearised rate maps ('lin') in case of linear track data. This function
            % is mainly to allow for flexibly choosing which parameters should be used
            % for the map construction.
            % For the parameter space check 'scanpix.maps.defaultParamsRateMaps'.
            %
            % Syntax:
            %       obj.addMaps
            %       obj.addMaps(mapType)
            %       obj.addMaps(mapType, trialInd)
            %       obj.addMaps(mapType, trialInd, varargin )
            %       obj.addMaps(mapType, [], varargin )
            %       obj.addMaps(__, 'load' )
            %       obj.addMaps(__, 'ui' )
            %       obj.addMaps(__, prmsStruct )
            %       obj.addMaps(__, Name-Value comma separated list )
            %
            % Inputs:
            %    mapType  - mapType string {'rate','dir', lin'}
            %    trialInd - numeric index for trial(s) for which you want to make maps (default: all trials)
            %    varargin - optional argument to control parameter selection
            %             - 'default': use default parameter space
            %             - 'load': load parameter file from disk
            %             - 'ui': open UI dialogue for parameter collection
            %             - prmsStruct: structure with parameter fields to be changed from defaults
            %             - name-value: comma separated list of name-value pairs
            %
            % Outputs:
            %
            % see also: scanpix.maps.makeRateMaps; scanpix.maps.makeDirMaps; scanpix.maps.makeLinRMaps; scanpix.maps.linearisePosData
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % quite a lot of input parsing/checking to do here %%%%%%%%
            if ~obj.loadFlag
                warning('scaNpix::ephys::addMaps: You need to load some data first before you can make any maps... Fairly obvious if you ask me.');
                return;
            end
            
            if nargin < 2
                str = {'rate','dir','lin','sac','objVect','speed'};
                [select, loadCheck] = listdlg('PromptString','Select what maps to make:','ListString',str,'ListSize',[160 100],'SelectionMode','Single');
                if ~loadCheck
                    warning('scaNpix::ephys::addMaps: No data selected. No maps will be created. Boring...');
                    return;
                else
                    mapType = str{select};
                end
            end
            
            if nargin < 3 || isempty(trialInd)
                trialInd = 1:length(obj.trialNames); % default: do all trials in object
            end
            
            if nargin < 4 || isempty(varargin{1})
                if isempty( obj.params('myRateMapParams') )
                    % use default params
                    prms = obj.mapParams.(mapType);
                else
                    % use user defined params
                    classFolder = what('scanpix');
                    try
                        temp    = load( fullfile(classFolder.path,'files', obj.params('myRateMapParams') ) ); %
                    catch
                        ME = MException('scaNpix::ephys::addMaps:''rateMapParamsNotFound', ['' fullfile(classFolder.path,'files', obj.params('myRateMapParams') ) ''' doesn''t exist buddy... ']);
                        throw(ME);
                    end
                    f            = fieldnames(temp);
                    prms         = temp.(f{1}).(mapType);
                    warning(['scaNpix::ephys::addMaps: Using user-defined parameters to make maps loaded from ' fullfile(classFolder.path,'files', obj.params('myRateMapParams') ) '.'] );
                end
            end
            
            if nargin == 4 && ischar(varargin{1})
                if strcmpi(varargin{1},'default')
                    % use defaults
                    prms = scanpix.maps.defaultParamsRateMaps;
                    prms = prms.(mapType);
                elseif strcmpi(varargin{1},'load')
                    % load from file
                    [fName, dataDir] = uigetfile(fullfile(classFolder.path,'files', '*.mat'), 'Select map params to load.');
                    % fail gracefully
                    if isnumeric(fName)
                        warning('scaNpix:ephys:addMaps: Loading Map Params aborted. Will use defaults instead!');
                        prms         = obj.mapParams.(mapType);
                    else
                        temp         = load( fullfile(dataDir, fName) );
                        f            = fieldnames(temp);
                        prms         = temp.(f{1}).(mapType);
                    end
                    
                elseif strcmpi(varargin{1},'ui')
                    % open UI dialogue to grab parameters
                    prms             = obj.mapParams.(mapType);
                    prompts          = fieldnames(prms);
                    defaultVals      = struct2cell(prms);
                    output           = scanpix.helpers.makeCustomUIDialogue(prompts, defaultVals);
                    % exit gracefully
                    if isempty(output)
                        warning('scaNpix::ephys::addMaps: Aborted changing parameters for rate map generation - will use those currently set in object.');
                        prms         = obj.mapParams.(mapType);
                    else
                        % format as structure
                        prms         = cell2struct(output(:,2), output(:,1), 1);
                    end
                    
                end
            end
            
            % Name-Value pair or struct input to change some params
            % explicitly
            if ~isempty(varargin) && (nargin > 4 || isstruct(varargin{1}) )
                prms = obj.mapParams.(mapType);
                if ischar(varargin{1})                                                           %
                    for i=1:2:length(varargin);   prms.(varargin{i}) = varargin{i+1};   end   %
                elseif isstruct(varargin{1})                                                     %
                    s = varargin{1};   f = fieldnames(s.(mapType));                                        %
                    for i=1:length(f);   prms.(f{i}) = s.(mapType).(f{i});   end                        %
                end
            end
            
            % these params get set from general data params container
            prms.posFs               = obj.params('posFs');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % make some maps
            switch lower(mapType)
                case 'rate'
                    prms.speedFilterLimits   = [prms.speedFilterLimitLow prms.speedFilterLimitHigh];
                    
                    for i = trialInd
                        if ~isempty( obj.trialMetaData(i).envSize )
                            prms.envSize = obj.trialMetaData(i).envSize / 100 * obj.trialMetaData(i).ppm; % in pixels
                        elseif strcmp(obj.fileType,'.set')
                            prms.envSize = [obj.trialMetaData(i).xmax-obj.trialMetaData(i).xmin obj.trialMetaData(i).ymax-obj.trialMetaData(i).ymin];
                        end
                        
                        [ obj.maps(1).rate{i}, obj.maps(1).pos{i}, obj.maps(1).spike{i} ] = scanpix.maps.makeRateMaps(obj.spikeData.spk_Times{i}, obj.posData.XY{i}, obj.spikeData.sampleT{i}, obj.trialMetaData(i).ppm, obj.posData.speed{i}, prms );
                    end
                    
                    %         % pad maps so size is the same for all maps is set
                    %         mapSz = reshape(cell2mat(cellfun(@(x) size(x{1}),obj.maps.rate,'uni',0)),2,[]);
                    %         maxSz = max(mapSz,[],2);
                    %         padSz = ceil((maxSz-mapSz) ./ 2);
                    %
                    %         mapType = {'rate','pos','spike'};
                    %         for i = 1:length(mapType)
                    %             for j = trialInd
                    %                 obj.maps.(mapType{i}){j} = cellfun(@(x) padarray(x,[padSz(1,j) padSz(2,j)],NaN), obj.maps.(mapType{i}){j},'uni',0);
                    %                 if any(mapSz(:,j) + 2*padSz(:,j) > maxSz)
                    %                     obj.maps.(mapType{i}){j} = cellfun(@(x) x(1:maxSz(1),1:maxSz(2)),obj.maps.(mapType{i}){j},'uni',0); % might need to trim off one row/column
                    %                 end
                    %             end
                    %         end
                    
                case 'dir'
                    prms.speedFilterLimits = [prms.speedFilterLimitLow prms.speedFilterLimitHigh];
                    for i = trialInd
                        obj.maps(1).dir{i} = scanpix.maps.makeDirMaps( obj.spikeData.spk_Times{i}, obj.posData.direction{i}, obj.spikeData.sampleT{i}, obj.posData.speed{i},  prms  );
                    end
                    
                case 'lin'
                    % map making a bit more involved mostly because we need some extra parameters we can't infer from data alone - should check if we
                    % can solve this more elegantly
                    prms.speedFilterLimits  = [prms.speedFilterLimitLow prms.speedFilterLimitHigh];
                    skipNextUI = false;
                    for i = trialInd
                        % we need the type of the track (for how smoothing is done)
                        if ~isempty(obj.trialMetaData(i).trialType)
                            trackProps.type = obj.trialMetaData(i).trialType;
                        else
                            uiInput = inputdlg({'linear track type', 'linear track length (cm)'},'',1,{'sqtrack','62.5'});
                            if isempty(uiInput)
                                warning('scaNpix::ephys::addMaps: No track properties, no linear rate maps...');
                                return;
                            else
                                trackProps.type = uiInput{1};
                                obj.trialMetaData(i).trialType = uiInput{1};
                                trackProps.length = str2double(uiInput{2});
                                obj.trialMetaData(i).trackLength = str2double(uiInput{2});
                                skipNextUI = true;
                            end
                        end
                        
                        % we need the length of the track for the pos scaling to work
                        if ~isempty(obj.trialMetaData(i).trackLength)
                            trackProps.length= obj.trialMetaData(i).trackLength;
                        elseif ~skipNextUI
                            uiInput = inputdlg({'linear track length (cm)'},'',1,{'62.5'});
                            if isempty(uiInput)
                                warning('scaNpix::ephys::addMaps: No track length, no linear rate maps...');
                                return;
                            else
                                trackProps.length = str2double(uiInput{1});
                                obj.trialMetaData(trialInd).trackLength = str2double(uiInput{1});
                            end
                        end
                        trackProps.ppm   = obj.trialMetaData(i).ppm;
                        trackProps.posFs = obj.params('posFs');
                        
                        [obj.maps(1).lin{i}, posMap, obj.posData(1).linXY{i}] = scanpix.maps.makeLinRMaps(obj.spikeData.spk_Times{i}, obj.posData.XY{i}, obj.spikeData.sampleT{i}, obj.posData.direction{i},obj.posData.speed{i}, trackProps, prms );
                        obj.maps(1).linPos{i} = num2cell(posMap,2);
                    end
                    
                case 'sac'
                    
                    if isempty(obj.maps(1).rate{1})
                        warning('scaNpix::ephys::addMaps: You need to generate rate maps before demanding spatial autocorrelograms.')
                        return
                    end
                    
                    for i = trialInd
                        % we don't want to smooth the AC but use smoothed rate map as input
                        obj.maps(1).sACs{i} = cellfun(@(x) scanpix.analysis.spatialCrosscorr(x,x),obj.maps.rate{i},'uni',0);
                    end
                    
                case 'objvect'
                    for i = trialInd
                        if isfield(obj.trialMetaData(i),'objectPos') && ~isempty(obj.trialMetaData(i).objectPos)
                            obj.maps(1).OV{i} = scanpix.maps.makeOVMap( obj.spikeData.spk_Times{i}, obj.posData.XY{i}, obj.spikeData.sampleT{i}, obj.trialMetaData(i).objectPos, obj.trialMetaData(i).ppm,  prms  );
                        else
                            warning('scaNpix::ephys::addMaps: If you to generate object vector maps, you need to have a field ''objectPos'' in your trialMetaData for any trial you want to generate these beauties for. Do that now and you won''t be disappointed');
                        end
                    end
                    
                case 'speed'
                    
                    for i = trialInd
                        obj.maps(1).speed{i} = scanpix.maps.makeSpeedMap( obj.spikeData.spk_Times{i}, obj.posData.speed{i}, obj.trialMetaData(i).duration,  prms  );
                    end
                    
                otherwise
                    ME = MException('scaNpix::ephys::addMaps:invalidMapType', ['' mapType ''' is not yet a supported map type. You need to invent that one yourself I am afraid... ']);
                    throw(ME);
            end
            
        end
        
        %%
        function spatProps = getSpatialProps(obj, score)
            
            switch lower(score)
                case {'spatialinfo','si'}
                    if isempty(obj.maps(1).rate) || isempty(obj.maps(1).rate{1})
                        warning('scaNpix::ephys::getSpatialProps::spatialInfo: You need to make rate maps before demanding their spatial info.');
                        spatProps = [];
                        return
                    end
                    
                    spatProps = nan(size(obj.cell_ID,1),length(obj.trialNames));
                    for i = 1:length(obj.trialNames)
                        spatProps(:,i) = scanpix.analysis.spatial_info(obj.maps(1).rate{i},obj.maps(1).pos{i});
                    end
                case {'raleighvect','rv'}
                    if isempty(obj.maps(1).dir) || isempty(obj.maps(1).dir{1})
                        warning('scaNpix::ephys::getSpatialProps::rVect:You need to make dir maps before demanding their rayleigh vector lengths.');
                        spatProps = [];
                        return
                    end
                    
                    spatProps = nan(size(obj.cell_ID,1),length(obj.trialNames));
                    for i = 1:length(obj.trialNames)
                        spatProps(:,i) = cell2mat(cellfun(@(x) scanpix.analysis.rayleighVect(x),obj.maps(1).dir{i},'uni',0));
                    end
                case {'gridprops','grid'}
                    if isempty(obj.maps(1).sACs) || isempty(obj.maps(1).sACs{1})
                        warning('scaNpix::ephys::getSpatialProps::gridProps:You need to make spatial ACs before demanding grid properties.');
                        spatProps = [];
                        return
                    end
                    
                    spatProps = nan(size(obj.cell_ID,1),5,length(obj.trialNames));
                    for i = 1:length(obj.trialNames)
                        [~, temp]        = cellfun(@(x) scanpix.analysis.gridprops(x,obj.mapParams.gridProps),obj.maps(1).sACs{i},'uni',0);
                        % for now just output the basics
                        spatProps(:,:,i) = cell2mat(cellfun(@(x) [x.gridness x.waveLength x.orientation],temp,'uni',0));
                    end
                case {'borderscore','bs'}
                    if isempty(obj.maps(1).rate) || isempty(obj.maps(1).rate{1})
                        warning('scaNpix::ephys::getSpatialProps::borderScore:You need to make rate maps before demanding their border score.');
                        spatProps = [];
                        return
                    end
                    
                    spatProps = nan(size(obj.cell_ID,1),length(obj.trialNames));
                    for i = 1:length(obj.trialNames)
                        [spatProps(:,i), ~] = scanpix.analysis.getBorderScore(obj.maps(1).rate{i},obj.mapParams.rate.binSizeSpat);
                    end
                otherwise
            end
        end
    end
    
    %%
    methods
        function loadWaves(obj,varargin)
            % loadWaves - load waveform data from neuropixel recordings
            % Method for scanpix class objects
            % Wrapper to load waveforms - in case you want to load every
            % single waveform, be prepared to wait a little...
            %
            % Syntax:  
            %    obj.loadWaves
            %    obj.loadWaves('ui')
            %    obj.loadWaves(Name-Value comma separated list )
            %
            % Inputs:
            %    varargin - 'ui' - bring up UI dialoge to select params
            %             - name-value: comma separated list of name-value
            %               pairs for params (see scanpix.npixUtils.extract_waveforms)
            %
            % Outputs:
            %
            % See also: scanpix.npixUtils.extract_waveforms
            %
            % LM 2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            if strcmp(obj.type,'dacq')
                warning('scaNpix::ephys::loadWaves: Waveforms for DACQ type objects are auto loaded. No need to ask for that again here compadre...');
                return
            end
            
            % deal with loading params
%             if nargin == 1
%                 % use defaults, only add channel n from object metadata
%                 addParams = {'nch'; obj.trialMetaData(1).nChan};
            if strcmp(varargin{1},'ui')
                % UI dialoge
%                 prompts = {'mode',    'N channels',                  'N waves / cluster', 'n chan / waveform', 'n samp / waveform', 'nSamplesPrePeak', 'apply CAR', 'save' };
%                 varargs = {'mode',    'nch',                         'nwave',              'getnch',           'nsamp',             'prepeak',         'car',       'save' };
%                 defVals = {'single',  obj.trialMetaData(1).nChanTot, 250,                  5,                  40,                  0.375,             0,           0,     };
                prompts = {'mode',    'N waves / cluster', 'n chan / waveform', 'n samp / waveform', 'nSamplesPrePeak', 'apply filter', 'save' };
                varargs = {'mode',    'nwave',              'getnch',           'nsamp',             'prepeak',         'filter',       'save' };
                defVals = {'single',  250,                  5,                  40,                  0.375,             0,              0,     };
                rtn = scanpix.helpers.makeCustomUIDialogue(prompts,defVals);
                if isempty(rtn)
                    warning('scaNpix::ephys::loadWaves: Waveform loading aborted. That lacks class mate...');
                    return;
                end

                addParams = cell(2,size(rtn,1));
                for i = 1:size(rtn,1)
                    addParams{1,i} = varargs{i};
                    addParams{2,i} = rtn{i,2};
                end
            elseif length(varargin) > 1
                % name-value pairs
                addParams = {varargin{1:2:end};varargin{2:2:end}};
%                 % always add channel n unless already supplied
%                 if ~any(strcmp(addParams(1,:),'nch'))
%                     addParams(:,end+1) = {'nch'; obj.trialMetaData(1).nChanTot};
%                 end
            end
            % extract waveforms
            scanpix.npixUtils.extract_waveforms(obj,1:length(obj.trialNames),addParams{:});
        end
        
    end
    
end



%             
%             if loadFromPhy    
%                 % load spike times
%                 spikeTimes         = readNPY(fullfile(obj.dataPath{trialIterator},'spike_times.npy'));
%                 spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
%                 obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
%                 % load cluster IDs
%                 clustIDs           = readNPY(fullfile(obj.dataPath{trialIterator},'spike_clusters.npy')) + 1; % 0 based index
%                 % phy curated output - only load clusters labelled good
%                 clu_info       = tdfread(fullfile(obj.dataPath{trialIterator},'cluster_info.tsv'),'tab');
%                 goodLabel      = strcmp(cellstr(clu_info.group),'good');
%                 try
%                     good_clusts    = clu_info.id(goodLabel) + 1;
%                 catch
%                     good_clusts    = clu_info.cluster_id(goodLabel) + 1;
%                 end
%                 clu_Depth      = clu_info.depth(goodLabel);  % I think this amd the following seem to not give correct results in somce cases?? Phy Issue?? 
%                 clu_Ch         = clu_info.ch(goodLabel) + 1; % this is 0 based 
%                 cluLabel       = string(clu_info.group);
%                 cluLabel       = cluLabel(goodLabel);
%             else
%                 if reloadFlag || exist(fullfile(obj.dataPath{trialIterator},'spike_times.npy'),'file') ~= 2
%                     % for a relaod we always request user interaction
%                     [fName, path2data_A] = uigetfile(fullfile(obj.dataPath{trialIterator},'*.npy'),['Select spike_times.npy for ' obj.trialNames{trialIterator}]);
%                     if isnumeric(fName)
%                         warning('Cannot continue reloading spike sorting results. Your journey ends here young padawan!');
%                         return
%                     end
%                     path2data_B = path2data_A;
%                 % try and load cluster IDs from back up folder if data was PHY'd as this will prevent issues further down in case you curated the data (and merged/split at least one cluster)    
%                 elseif exist(fullfile(obj.dataPath{trialIterator},'cluster_info.tsv'),'file') == 2
%                     
%                     if isfolder(fullfile(obj.dataPath{trialIterator},'backUpFiles'))
%                         path2data_A = fullfile(obj.dataPath{trialIterator},'backUpFiles');
%                         path2data_B = obj.dataPath{trialIterator};
%                     else
%                         [path2data_A,path2data_B] = deal(obj.dataPath{trialIterator});
%                         warning('scaNpix::npix::loadSpikes:Can''t find folder with BackUp files! If you merged/split clusters in PHY, loading raw kilosort results will fail shortly...');
%                     end
%                 else
%                     [path2data_A,path2data_B] = deal(obj.dataPath{trialIterator});
%                 end
                
                
                % load spike times
%                 spikeTimes         = readNPY(fullfile(path2data_A,'spike_times.npy'));
%                 spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
%                 obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
%                 %
%                 clustIDs  = readNPY(fullfile(path2data_A,'spike_clusters.npy')) + 1; % 0 based index
%                 ks_labels = tdfread(fullfile(path2data_A,'cluster_KSLabel.tsv'),'tab'); % we need the raw kilosort output here before you ran Phy as this file gets overwritten! - maybe point to backup folder?
%                 
%                 try
%                     % try and load cluster IDs from back up folder as this will prevent issues further down in case you curated the data in phy (and merged/split at least one cluster)
%                     % load spike times
%                     spikeTimes         = readNPY(fullfile(obj.dataPath{trialIterator},'backUpFiles','spike_times.npy'));
%                     spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
%                     obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
%                     %
%                     clustIDs  = readNPY(fullfile(obj.dataPath{trialIterator},'backUpFiles','spike_clusters.npy')) + 1; % 0 based index
%                     ks_labels = tdfread(fullfile(obj.dataPath{trialIterator},'backUpFiles','cluster_KSLabel.tsv'),'tab'); % we need the raw kilosort output here before you ran Phy as this file gets overwritten! - maybe point to backup folder?
%                 catch
%                     warning('scaNpix::npix::loadSpikes:Can''t find folder with BackUp files! If you PHY''d the data and merged/split clusters, loading raw kilosort results will fail shortly...');
%                     % load spike times
%                     spikeTimes         = readNPY(fullfile(obj.dataPath{trialIterator},'spike_times.npy'));
%                     spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
%                     obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
%                     % load cluster IDs
%                     clustIDs  = readNPY(fullfile(obj.dataPath{trialIterator},'spike_clusters.npy')) + 1; % 0 based index
%                     ks_labels = tdfread(fullfile(obj.dataPath{trialIterator},'cluster_KSLabel.tsv'),'tab'); % we need the raw kilosort output here before you ran Phy as this file gets overwritten! - maybe point to backup folder?
%                 end
                % raw kilosort output - we'll take everything in this case and can
                % filter later if needed
%                 cluLabel       = string(ks_labels.KSLabel); % this is either 'good' or 'mua'
%                 %     ind_good       = strcmp(cluLabel,'good') | strcmp(cluLabel,'mua');
%                 good_clusts    = ks_labels.cluster_id + 1; % ID for clusters also 0 based
%                 % get depth estimate for each cluster based on template ampl
%                 % distribution across probe
%                 templates          = readNPY(fullfile(path2data_B,'templates.npy'));
%                 Winv               = readNPY(fullfile(path2data_B,'whitening_mat_inv.npy'));
%                 chanPos            = readNPY(fullfile(path2data_B,'channel_positions.npy'));
%                 chanMapKS          = double(readNPY(fullfile(path2data_B,'channel_map.npy'))) + 1;  % 
%                 [clu_Depth,clu_Ch] = scanpix.npixUtils.getCluChDepthFromTemplates(templates, Winv, [chanMapKS(:) chanPos(:,2)]);
%             end
