classdef nexus < handle
    % nexus data class. Let's you load in neuronexus data into a class
    % object for analysis or data inspection
    %
    % Class constructor: scanpix.nexus
    % Construct class obj
    %
    % Usage:
    %       obj = scanpix.nexus;
    %       obj = scanpix.nexus(prmsMode);
    %       obj = scanpix.nexus(prmsMode, uiFlag);
    %
    % Inputs:
    %       prmsMode: 'default' - uses default parameter (default)
    %                 'ui'      - opens UI dialogue
    %                 'file'    - load params from file
    %
    %         uiFlag: true (default) or false
    %                 - if false will skip UI dialogue for data selection
    %                   (e.g. when you use constructor programmatically)
    %
    %  IV
    %% PROPERTIES
    properties % params % -from npix obj
        % Map object
        params                containers.Map
        chanMap               struct
    end
    
    properties % meta data % -from npix obj
        dataPath(1,:)         string
        dataSetName(1,:)      char
        trialNames(1,:)       string
        cell_ID(:,4)          double %{mustBePositive, mustBeNonNan, mustBeNonzero}
        cell_Label(:,1)       string
    end
    
    properties % trial data % -from dacq obj
        trialMetaData(1,:)    struct  % leave non-scalar as many fields and indexing into it, is not required very often
        posData               struct  = struct('XYraw',[],'XY',[],'direction',[],'speed',[],'linXY',[]);
        spikeData             struct  = struct('spk_Times',[],'spk_waveforms',[],'sampleT',[]);
        lfpData               struct  = struct('lfp',[],'lfpHighSamp',[],'lfpTet',[]);% this is from dacq - need to decide where I'm getting lfp from - IV
    end

    properties % maps % - from dacq obj
        maps                  struct  = struct('rate',[],'spike',[],'pos',[],'dir',[],'sACs',[],'OV',[],'speed',[]);
        linMaps               struct  = struct('linRate',[],'linPos',[],'linRateNormed',[]);
    end
    
    properties(Hidden)
        fileType              char    = '.set';
        uniqueCellIndex(:,1)  logical % this line is commented out in the npix bit 
        fields2spare          cell    = {'params','dataSetName','cell_ID','cell_Label'};  % spare this when deleting or rearranging data. make sure to add new properties that should be spared here!
        mapParams             struct  = scanpix.maps.defaultParamsRateMaps; 
        loadFlag              logical = false;                               % flag so we know something has been loaded into object
        badChans              %npix only 
    end
    %% METHODS
    
    % constructor %
    methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = nexus(prmsMode,uiFlag)
            % nexus - creates class object 
            %
            % Syntax:  
            %       obj = nexus;
            %       obj = nexus(prmsMode);
            %       obj = nexus(prmsMode, uiFlag);
            %
            % Inputs:
            %    prmsMode - 'default' - uses default parameter (default)
            %               'ui'      - opens UI dialogue
            %               'file'    - load params from file
            %    uiFlag   - true (default)/false - skip UI set file
            %               selection dialogue if false
            %
            % Outputs:
            %    obj      - nexus object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            if nargin == 0
                obj.params = scanpix.helpers.getParams(obj,'default');
                uiFlag     = true;
            else
                obj.params = scanpix.helpers.getParams(obj, prmsMode);
                if nargin < 2
                    uiFlag     = true;
                end
            end
            
            if isempty(obj.params)
                warning('scaNpix: No params selected. Using defaults.');
                obj.params = scanpix.helpers.defaultParamsContainer(class(obj));
            end

            if uiFlag
                [obj.trialNames, obj.dataPath] = scanpix.helpers.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
                if isempty(obj.trialNames); return; end
                scanpix.helpers.selectTrials(obj.trialNames, obj);
            end
            
            if ~isempty(obj.dataPath)
                obj.params('defaultDir') = [fileparts(obj.dataPath{1}) filesep];  % use parent directory rather than absolut as this makes more sense for multiple trial data
            end
        end
    end
     %% took this whole section from npix
    
    % data loading %
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
            % see also: obj.loadMeta; obj.loadPos; obj.loadSpikes; obj.loadLFPs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isempty(obj.dataPath)
                [obj.trialNames, obj.dataPath] = scanpix.helpers.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
                scanpix.helpers.selectTrials(obj.trialNames, obj);
            end
            
            if nargin < 2
                str = {'all','pos','spikes','lfp'};
                [select, loadCheck] = listdlg('PromptString','Select what data to load:','ListString',str,'ListSize',[160 100]);
                if ~loadCheck
                    warning('scaNpix::load: No data selected. Nothing is loaded. Boring...');
                    return;
                else
                    loadMode = str(select);
                end
            end
            
            if ~iscell(loadMode)
                loadMode  = {loadMode};
            end
            
            if nargin < 3
                loadStr = obj.trialNames;
            else
                loadStr = obj.trialNames( ismember(obj.trialNames,varargin{1}) );
            end
            
            trialInd = find(ismember(obj.trialNames,loadStr));
            
            for i = 1:length(loadStr)
                
                obj.loadMeta(trialInd(i)); % load some meta data - will need to make a meta data file compatible with LM's code -IV
                               
                for j = 1:length(loadMode)
                    
                    switch lower(loadMode{j})
                        case 'all'
                            obj.loadPos(trialInd(i));
                            obj.loadSpikes(trialInd(i));
                            %                             obj.loadLFPs(trialInd(i));
                            break;
                        case 'pos'
                            obj.loadPos(trialInd(i));
                        case 'spikes'
                            obj.loadSpikes(trialInd(i));
                        case 'lfp'
                            %%% TO DO !!!!!
                            %                             obj.loadLFPs(trialInd(i));
                        case 'meta'
                            % already loaded!
                        otherwise
                            ME = MException('scaNpix:load:invalidInput', ['' loadMode{j} ''' is not a valid data type selector. Try ''all'', ''pos'', ''spikes'' and/or ''eeg'' instead. ']);
                            throw(ME);
                    end
                end
            end
            
            %             obj.metaData = [];
            %             obj.addMetaData(obj.defaultMetaDataFields); % add the default meta data fields  %%% XML FILE SHOULD CONTAIN THIS!
            %
            obj.loadFlag = true; % flag loading
            
            % pre-allocate some data (e.g. part that isn't loaded), such
            % that all properties have same size
            scanpix.helpers.preallocEmpty(obj,true,{'posData','spikeData','lfpData'});
        end
        
    end
end 