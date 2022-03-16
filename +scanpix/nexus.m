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
     %% dacq loading functions
    methods(Hidden)
        %%
        function loadSet(obj, trialIterator)
            % loadSet - load set file data from DACQ files
            % Method for dacq class objects (hidden)
            %
            % Syntax:  obj.loadSet(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %
            % See also:
            %
            % TW/LM 2020 (adapted from org. SCAN function)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % sanity check
            if ~exist(fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.set']),'file')
                ME = MException('scaNpix:loadSet:setFileNotFound', ['Could not find ''' strrep(obj.dataPath{trialIterator},'\','/') obj.trialNames{trialIterator} '.set''' ]); % should maybe switch to forward slash everywhere
                throw(ME);
            end
            
            % Read file %
            fid = fopen(fullfile(obj.dataPath{trialIterator}, [obj.trialNames{trialIterator} '.set']), 'r');
            C = textscan(fid, '%s %[^\r\n]');
            sFileTxt = [cat(1,C{1}) cat(1,C{2})];
            fclose(fid);
            
            % legacy data isn't supported. multiplex pre-amp has different%
            % ADC range - change this for multiplex!! -IV
            if str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ADC_fullscale_mv')) ~= 1500 && isempty(scanpix.dacqUtils.getValue(sFileTxt, 'demux_en_dac_1')) 
                ME = MException('scaNpix:loadSet:legacyDataError', ['''' strrep(obj.dataPath{trialIterator},'\','/') obj.trialNames{trialIterator} '.set'' was not acquired with DACQ USB. Data obtained from legacy versions are not supported.' ]);
                throw(ME);
            end
            
            %%% These are straightforward one value fields %%%
            setFile.tracked_spots    = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'tracked_spots'));
            setFile.xmin             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'xmin'));
            setFile.xmax             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'xmax'));
            setFile.ymin             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ymin'));
            setFile.ymax             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ymax'));
            setFile.sw_version       = scanpix.dacqUtils.getValue(sFileTxt, 'sw_version');   % Don't STR2DOUBLE - can be 4.00a
            setFile.trial_time       = scanpix.dacqUtils.getValue(sFileTxt, 'trial_time');
            setFile.ADC_fullscale_mv = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ADC_fullscale_mv'));
            duration                 = sscanf(scanpix.dacqUtils.getValue(sFileTxt, 'duration'),'%d');
            if ~~rem(duration,10)
                duration             = duration - rem(duration,10);
            end
            setFile.duration         = duration;
            
            % light parameters: make a 1:4 vector %
            lightParams = {'lightBearing', 'colactive'};
            for i = 1:2
                for j = 1:4
                    setFile.(lightParams{i})(j) = str2double(scanpix.dacqUtils.getValue(sFileTxt, [lightParams{i} '_' num2str(j)]));
                end
            end
            
            % Gains %
            for i = 1:128 % HARD CODED (but Jim will probably never make a higher channel count system...)
                setFile.gains(i) = str2double(scanpix.dacqUtils.getValue(sFileTxt, ['gain_ch_' num2str(i-1)]));
            end
            setFile.gains     = reshape(setFile.gains,4,32)';
            setFile.fullscale = (setFile.ADC_fullscale_mv ./ setFile.gains) .* 1000;
            
            %%% EEG Channels %%%
            % Which channels EEGs recorded? %
            recordingChannel = zeros([1 128]); %HARD CODED
            for i = 1:length(recordingChannel) %%
                temp=scanpix.dacqUtils.getValue(sFileTxt,['saveEEG_ch_' num2str(i)]);
                if isempty(temp)
                    break;
                elseif str2double(temp) % temp is '1' or '0' for EEG used or not.
                    recordingChannel(i) = str2double(scanpix.dacqUtils.getValue(sFileTxt,['EEG_ch_'  num2str(i)]));
                end
            end
            EEGSlotActive    = find(recordingChannel);                  % Now a list of eeg channels in use ..
            recordingChannel = recordingChannel(recordingChannel~=0);   %  .. and their recording slots. (recordingChannel=0 if slot not in use).
            
            % Get signal source channel, gain, filters %
            if isempty(recordingChannel)
                [setFile.lfp_channel,setFile.lfp_recordingChannel,setFile.lfp_slot,setFile.lfp_scalemax,setFile.lfp_filter,setFile.lfp_filtresp,setFile.lfp_filtkind,setFile.lfp_filtfreq1,...
                    setFile.filtfreq2,setFile.lfp_filtripple] = deal([]); % In case of null EEG
            else
                for i=1:length(recordingChannel)
                    mode = str2double(scanpix.dacqUtils.getValue(sFileTxt, ['mode_ch_' num2str(recordingChannel(i)-1)]));
                    if any( mode == [1 3] ) % Mode B or -B (ref or -ref)
                        chTemp = scanpix.dacqUtils.getValue(sFileTxt, ['b_in_ch_' num2str(recordingChannel(i)-1)]);
                        % if not free referencing system need to fetch ref channel ID from ref field
                        if ~str2double(scanpix.dacqUtils.getValue(sFileTxt,'modeanalog32'))
                            chTemp = scanpix.dacqUtils.getValue(sFileTxt, ['ref_' chTemp]);
                        end
                        setFile.lfp_channel(i) = str2double(chTemp) + 1;
                    else
                        setFile.lfp_channel(i) = recordingChannel(i);
                    end
                    setFile.lfp_recordingChannel(i) = recordingChannel(i);   % Actual recording channel, not signal source (useful when approaching raw data).
                    setFile.lfp_slot(i)             = EEGSlotActive(i);      % Is this EEG1, EEG2, etc (according to the 'record>setup>EEG' tab.
                    setFile.lfp_scalemax(i)         = setFile.fullscale( ceil(recordingChannel(i)/4), recordingChannel(i)-(((ceil(recordingChannel(i)/4))-1)*4) );
                    setFile.lfp_filter(i)           = str2double(scanpix.dacqUtils.getValue(sFileTxt, ['filter_ch_' num2str(recordingChannel(i)-1)]));
                    
                    % DacqUSB only - fuller filter spec %
                    fSpec = {'filtresp', 'filtkind', 'filtfreq1', 'filtfreq2', 'filtripple'};
                    for j = 1:length(fSpec)
                        setFile.(['lfp_' fSpec{j}]) = str2double(scanpix.dacqUtils.getValue(sFileTxt, [fSpec{j} '_ch_' num2str(recordingChannel(i)-1)]));
                    end
                    
                end
            end
            % need to initialise here
            setFile.ppm         = []; % will be populated when loading pos data
            setFile.ppm_org     = []; % will be populated when loading pos data
            setFile.trialType   = []; % no method yet for dacq data to gather that - integrate metadata maybe a .csv file at this stage IV 
            setFile.trackLength = []; % no method yet for dacq data to gather that
            setFile.envSize     = []; % no method yet for dacq data to gather that
            
            % output
            if isempty(obj.trialMetaData)
                obj.trialMetaData = setFile;
            else
                obj.trialMetaData(trialIterator) = setFile;
            end

            % load channel map -this bit is from npix 
            chanMapInfo = dir(fullfile(obj.dataPath{trialIterator},'*kiloSortChanMap*'));
            if isempty(chanMapInfo)
                chanMapName = scanpix.npixUtils.SGLXMetaToCoords_v2(fullfile(obj.dataPath{trialIterator})); % create map on the fly
            else
                chanMapName = fullfile(obj.dataPath{trialIterator}, chanMapInfo.name);
            end
            chanMapStruct = load(chanMapName);
            f = fieldnames(chanMapStruct);
            for i = 1:length(f)
                obj.chanMap(trialIterator).(f{i})  = chanMapStruct.(f{i});
            end
            obj.trialMetaData(trialIterator).nChan = sum(chanMapStruct.connected);
            
            %back to dacq from here
            
            if isempty(obj.dataSetName) %might make sense to use the npix one for this bit -IV keep track of what's happening
                obj.dataSetName = ['r_' obj.trialNames{trialIterator}(1:6)]; %%%%% THIS NEEDS FIX ONCE WE ARE CLEAR HOW TO HANDLE METADATA IN DACQ - make some sort of .csv containing necessary trial metadata - IV
            end
            
        end 
            %%
        function loadPos(obj, trialIterator)
            % loadPos - load position data from DACQ files
            % Method for dacq class objects (hidden)
            %
            % Syntax:  obj.loadPos(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %
            % See also:
            %
            % TW/LM 2020 (adapted from org. SCAN function)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('Loading pos file for %s [..........] ', obj.trialNames{trialIterator});
            
            %%% Read data from File %%%
            fid = fopen(fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.pos']),'r','ieee-be');
            if fid ==-1
                ME = MException('scaNpix::load_pos:posFileNotFound', ['Could not open pos file ''' fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.pos']) '']);
                throw(ME);
            end
            
            C = textscan(fid,'%s %[^\r\n]', 27);  % CAUTION! N Lines of header hard-coded.
            frewind(fid);
            posHeader    = [cat(1,C{1}) cat(1,C{2})];
            % Look for data start marker
            headerText   = fread(fid,800,'int8');
            headerText   = char(abs(headerText))';
            ds           = strfind(headerText,'data_start') + 10;
            % Read data %
            fseek(fid,(ds-1),'bof');                                    % Set file position indicator at data start
            nPosSamples  = sscanf(scanpix.dacqUtils.getValue(posHeader,'num_pos_samples'), '%d');
            data         = fread(fid, nPosSamples * 10, 'uint16');
            fclose(fid);
            
            %%% Reshape into correct output format %%%
            data         = reshape(data, [10 nPosSamples])';
            data         = data(:,3:10);
            data         = reshape(data, [nPosSamples, 2, 4]); % Now in format [nSamp, x:y, nLight]
            % Separate numpix, if existing, switch format to [nSamp, nLight, x:y] %
            nLights      = sum(obj.trialMetaData(trialIterator).colactive);
            if strcmp(scanpix.dacqUtils.getValue(posHeader,'pos_format'), 't,x1,y1,x2,y2,numpix1,numpix2') && nLights <= 2
                led_pos  = nan(nPosSamples, nLights, 2);
                led_pix  = nan(nPosSamples, nLights);
                for i = 1:nLights
                    for j = 1:2
                        led_pos(:,i,j) = data(:,j,i);
                    end
                    led_pix(:,i) = data(:,i,3); % numpix always seems to start at 3rd light (5th entry)
                end
            else
                % not sure this is still necessary as we are not allowing legacy data
                % anyway
                led_pos = nan(nPosSamples, nLights, 2);
                led_pix = [];
                for i = 1:nLights
                    for j = 1:2
                        led_pos(:,i,j) = data(:,j,i);
                    end
                end
            end
            led_pos(led_pos==1023) = NaN;   % mTint functions are
            led_pix(led_pos==1023) = NaN;   % expecting this.
            
            % scale pos data to specific PPM
%             obj.params('ppm_org') = sscanf(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'),'%d');
            obj.trialMetaData(trialIterator).ppm_org = sscanf(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'),'%d');

            if ~isempty(obj.params('ScalePos2PPM'))
                scaleFact = (obj.params('ScalePos2PPM') / str2double(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre')));
                led_pos = floor( led_pos .*  scaleFact);
                posHeader{strcmp(posHeader(:,1),'pixels_per_metre'),2} = num2str(obj.params('ScalePos2PPM'));
                obj.trialMetaData(trialIterator).xmin = obj.trialMetaData(trialIterator).xmin * scaleFact;
                obj.trialMetaData(trialIterator).xmax = obj.trialMetaData(trialIterator).xmax * scaleFact;
                obj.trialMetaData(trialIterator).ymin = obj.trialMetaData(trialIterator).ymin * scaleFact;
                obj.trialMetaData(trialIterator).ymax = obj.trialMetaData(trialIterator).ymax * scaleFact;
            end
%             obj.params('ppm') = sscanf(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'),'%d');
            obj.trialMetaData(trialIterator).ppm = sscanf(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'),'%d');

            % post process
            posFile.header  = posHeader; % For mTint function compatibility.
            posFile.led_pos = led_pos;
            posFile.led_pix = led_pix;
            %Interpolate, Filter and Smooth %  %%%%%%%%%% THIS COULD USE A BIT
            %MORE EFFICIENCY AND BE UPDATED
            [obj.posData(1).XYraw{trialIterator}, obj.posData(1).direction{trialIterator}, obj.posData(1).speed{trialIterator}] = scanpix.dacqUtils.postprocess_pos_data_v2(posFile, obj.params('posMaxSpeed'), obj.params('posSmooth'), obj.trialMetaData(trialIterator), obj.params('posHead'));
            
            obj.params('posFs') = sscanf(scanpix.dacqUtils.getValue(posHeader,'sample_rate'),'%d');
            % remove DACQ overhang
            if obj.trialMetaData(trialIterator).duration * obj.params('posFs') < length(led_pos)
                obj.posData(1).XYraw{trialIterator}     = obj.posData(1).XYraw{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs'),:); % truncate data
                obj.posData(1).direction{trialIterator} = obj.posData(1).direction{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs')); % truncate data
                obj.posData(1).speed{trialIterator}     = obj.posData(1).speed{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs')); % truncate data
            end
            
            % convert to integers
            obj.posData(1).XY{trialIterator} = [double( floor(obj.posData(1).XYraw{trialIterator}(:,1)) + 1 ), double( floor(obj.posData(1).XYraw{trialIterator}(:,2)) + 1 )];  %%% NECESSARY??
            
            
            fprintf('  DONE!\n');
            
        end
     %% - this bit is from npix - need to compatibilize - IV
        function loadSpikes(obj, trialIterator)
            % loadSpikes - load spike data from neuropixel files
            % Method for scanpix class objects (hidden)
            % We will just load spike times
            %
            % Syntax:  obj.loadSpikes(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %
            % See also: scanpix;
            %
            % LM 2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% To-Do
            
            % initialise
            fprintf('Loading Neural Data for %s .......... ', obj.trialNames{trialIterator});
            
            %% grab sync channel data
            cd(obj.dataPath{trialIterator});
            if ~isfield(obj.trialMetaData,'BonsaiCorruptFlag')
                syncTTLs = scanpix.npixUtils.loadSyncData();
            else
                syncTTLs = scanpix.npixUtils.loadSyncData(length(obj.posData.sampleT{trialIterator}),obj.trialMetaData(trialIterator).BonsaiCorruptFlag); 
            end
           
            
            % decide what to load - phy or kilosort          
            if obj.params('loadFromPhy')
                if exist(fullfile(obj.dataPath{trialIterator},'cluster_info.tsv'),'file') == 2
                    loadFromPhy = true;
                else
                    warning('scaNpix::npix::loadSpikes:Can''t find ''cluster_info'' from phy output. Will try using kilosort data instead!');
                    loadFromPhy = false;
                end
            else
                loadFromPhy = false;
            end
            
            
            %% load sorting resuts
            
            if loadFromPhy    
                % load spike times
                spikeTimes         = readNPY(fullfile(obj.dataPath{trialIterator},'spike_times.npy'));
                spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
                obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
                % load cluster IDs
                clustIDs           = readNPY(fullfile(obj.dataPath{trialIterator},'spike_clusters.npy')) + 1; % 0 based index
                % phy curated output - only load clusters labelled good
                clu_info       = tdfread(fullfile(obj.dataPath{trialIterator},'cluster_info.tsv'),'tab');
                goodLabel      = strcmp(cellstr(clu_info.group),'good');
                try
                    good_clusts    = clu_info.id(goodLabel) + 1;
                catch
                    good_clusts    = clu_info.cluster_id(goodLabel) + 1;
                end
                clu_Depth      = clu_info.depth(goodLabel);  % I think this amd the following seem to not give correct results in somce cases?? Phy Issue?? 
                clu_Ch         = clu_info.ch(goodLabel) + 1; % this is 0 based 
                cluLabel       = string(clu_info.group);
                cluLabel       = cluLabel(goodLabel);
            else
                try
                    % try and load cluster IDs from back up folder as this will prevent issues further down in case you curated the data in phy (and merged/split at least one cluster)
                    % load spike times
                    spikeTimes         = readNPY(fullfile(obj.dataPath{trialIterator},'backUpFiles','spike_times.npy'));
                    spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
                    obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
                    %
                    clustIDs  = readNPY(fullfile(obj.dataPath{trialIterator},'backUpFiles','spike_clusters.npy')) + 1; % 0 based index
                    ks_labels = tdfread(fullfile(obj.dataPath{trialIterator},'backUpFiles','cluster_KSLabel.tsv'),'tab'); % we need the raw kilosort output here before you ran Phy as this file gets overwritten! - maybe point to backup folder?
                catch
                    warning('scaNpix::npix::loadSpikes:Can''t find folder with BackUp files! If you PHY''d the data and merged/split clusters, loading raw kilosort results will fail shortly...');
                    % load spike times
                    spikeTimes         = readNPY(fullfile(obj.dataPath{trialIterator},'spike_times.npy'));
                    spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
                    obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
                    % load cluster IDs
                    clustIDs  = readNPY(fullfile(obj.dataPath{trialIterator},'spike_clusters.npy')) + 1; % 0 based index
                    ks_labels = tdfread(fullfile(obj.dataPath{trialIterator},'cluster_KSLabel.tsv'),'tab'); % we need the raw kilosort output here before you ran Phy as this file gets overwritten! - maybe point to backup folder?
                end
                % raw kilosort output - we'll take everything in this case and can
                % filter later if needed
                cluLabel       = string(ks_labels.KSLabel); % this is either 'good' or 'mua'
                %     ind_good       = strcmp(cluLabel,'good') | strcmp(cluLabel,'mua');
                good_clusts    = ks_labels.cluster_id + 1; % ID for clusters also 0 based
                % get depth estimate for each cluster based on template ampl
                % distribution across probe
                templates          = readNPY(fullfile(obj.dataPath{trialIterator},'templates.npy'));
                Winv               = readNPY(fullfile(obj.dataPath{trialIterator},'whitening_mat_inv.npy'));
                chanPos            = readNPY(fullfile(obj.dataPath{trialIterator},'channel_positions.npy'));
                chanMapKS          = double(readNPY(fullfile(obj.dataPath{trialIterator},'channel_map.npy'))) + 1;  % 
                [clu_Depth,clu_Ch] = scanpix.npixUtils.getCluChDepthFromTemplates(templates, Winv, [chanMapKS(:) chanPos(:,2)]);
            end
            
            % now we need to remove bad clusters and spike times outside trial
            unClustIDs         = unique(clustIDs);
            % index for 'good' clusters
            unGoodClustIDs = unClustIDs(ismember(unClustIDs,good_clusts)); % remove 'mua' or 'noise' clusters from list in case we deal with phy output
            [~,indGood]    = ismember(clustIDs,unGoodClustIDs); % only keep these
            % index for spike times within ]1st TTL, trialduration]
            indSTimes      = spikeTimes > 0 & spikeTimes <= obj.trialMetaData(trialIterator).duration;
            % final index of spike times we want to keep
            ind2keep       = indGood ~= 0 & indSTimes;
            
            % remove from data
            clustIDs       = clustIDs(ind2keep);
            spikeTimes     = spikeTimes(ind2keep);
            
            % sort clusters, so accumarray output is sorted
            [clustIDs, sortInd] = sort(clustIDs);
            spikeTimes          = spikeTimes(sortInd);
            
            % reformat into more convenient form
            spikeTimesFin  = accumarray(clustIDs,spikeTimes,[max(unGoodClustIDs) 1],@(x) {x});
            % remove empty clusters - need to make sure not to remove cells that
            % only fire in some trials of a trial sequence (we are assuming here that
            % you clustered all data together and then split back into individual
            % trials)
            OtherClusters    = true(length(spikeTimesFin),1);
            OtherClusters(good_clusts) = false;
            indEmpty       = cellfun('isempty',spikeTimesFin) & OtherClusters;
            spikeTimesFin  = spikeTimesFin(~indEmpty);
            
            %%%%%%%% I SHOULD CHECK WHETHER THIS WORKS FOR ALL CIRCUMSTANCES!!!!!!!
            if ~loadFromPhy
                clu_Ch         = clu_Ch(~indEmpty);
                clu_Depth      = clu_Depth(~indEmpty);
            end
            % sort by depth
            [clu_Depth, indSort] = sort(clu_Depth,'ascend');
            spikeTimesFin        = spikeTimesFin(indSort);
            %% output
            obj.spikeData(1).spk_Times{trialIterator} = spikeTimesFin;
            endIdxNPix                                = min( [ length(syncTTLs), find(syncTTLs < obj.trialMetaData(trialIterator).duration + syncTTLs(1),1,'last') + 1]);
            obj.spikeData(1).sampleT{trialIterator}   = syncTTLs(1:endIdxNPix) - syncTTLs(1);
            % we'll only need to grab this once - this assumes data was clustered together, but wouldn't make much sense otherwise?
            % what about reload?
            if trialIterator == 1
                good_clusts          = good_clusts(indSort);
                cluLabel             = cluLabel(indSort);
                clu_Ch               = clu_Ch(indSort);
                % likely at least the ref channel will have been removed before sorting - this will map channel ID back to raw data
                clu_Ch_mapped = scanpix.npixUtils.mapChans(obj.chanMap.connected,clu_Ch);

                obj.cell_ID    = [good_clusts, clu_Depth, clu_Ch clu_Ch_mapped];
                obj.cell_Label = cluLabel;
            end
            
            % also keep only tracking data within trial times
            if ~isempty(obj.posData.XYraw) && ~isempty(obj.posData.XYraw{trialIterator})
                obj.posData.XYraw{trialIterator}     = obj.posData.XYraw{trialIterator}(1:endIdxNPix,:);
                obj.posData.XY{trialIterator}        = obj.posData.XY{trialIterator}(1:endIdxNPix,:);
                obj.posData.speed{trialIterator}     = obj.posData.speed{trialIterator}(1:endIdxNPix,:);
                obj.posData.direction{trialIterator} = obj.posData.direction{trialIterator}(1:endIdxNPix,:);
                obj.posData.sampleT{trialIterator}   = obj.posData.sampleT{trialIterator}(1:endIdxNPix,:);
            end
            
            fprintf('  DONE!\n');
        end
    %% this bit is from dacq 
        function loadLFPs(obj, trialIterator)
            % loadLFPs - load eeg data from DACQ files
            % Method for dacq class objects (hidden)
            % We will load all low and (if available) high sample rate EEGs and we also fetch the tetrode ID each eeg was recorded from
            %
            % Syntax:  obj.loadLFPs(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %
            % See also: dacq;
            %
            % TW/LM 2020 (adapted from org. SCAN function)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            clear scanpix.fxchange.textprogressbar  % get rid of persistent variable in 'textprogressbar' which can cause issues if you stop code during execution
            
            noEGFflag = false;
            lfp2load = scanpix.dacqUtils.findDACQFiles([obj.dataPath{trialIterator} filesep], obj.trialNames{trialIterator},'eeg');
            if ischar(lfp2load)
                lfp2load = {lfp2load}; % in case only 1 eeg
            end
            
            sRateStr = {'lfpFs','lfpHighFs'};
            
            for i = 1:2
                % switch extension
                if i == 2
                    lfp2load = strrep(lfp2load,'.eeg','.egf');
                end
                
                scanpix.fxchange.textprogressbar(['Loading ' lfp2load{1}(end-3:end) ' data for ' obj.trialNames{trialIterator} ' ']);
                
                for j = 1:length(lfp2load)
                    
                    fid = fopen(fullfile(obj.dataPath{trialIterator}, lfp2load{j}),'r','ieee-be');  % 'ieee-be' is machine format, 'big endian'.
                    if fid == -1
                        noEGFflag = true;
                        scanpix.fxchange.textprogressbar(0);
                        continue
                    end
                    
                    hdr = fread(fid,400,'int8');
                    ds  = strfind(hdr','data_start') + 10; % data start marker
                    % more convenient format - read some info from header
                    frewind(fid);
                    %%%%%%%
                    tempHeader = textscan(fid,'%s %[^\n]',11);
                    tempHeader = horzcat(tempHeader{:});
                    
                    %%%%%%%
                    sRateInd = strcmp('sample_rate',tempHeader(:,1));
                    if sum( sRateInd ) == 0   % Sometimes there are 'empty' EEG files, which don't even have a full header. If we have one of these, just quit now.
                        warning(['scaNpix: Problem loading EEG. Empty EEG data file for ' lfp2load{j}]);
                        continue
                    end
                    obj.params(sRateStr{i}) = sscanf(tempHeader{sRateInd,2},'%d');  % hard code??
                    %%%%%%
                    bytesInd     = strcmp('bytes_per_sample',tempHeader(:,1));
                    bytesPerSamp = sscanf(tempHeader{bytesInd,2},'%d');
                    % read data
                    fseek(fid,ds,'bof');
                    if bytesPerSamp == 1
                        nSamplesInd    = strcmp('num_EEG_samples',tempHeader(:,1));
                        nSamples       = sscanf(tempHeader{nSamplesInd,2},'%d');
                        tempData       = fread(fid,nSamples,'int8');
                        obj.lfpData(1).lfp{trialIterator}{j}  = (double(tempData)./2^7) .* obj.trialMetaData(trialIterator).lfp_scalemax(j); %voltages
                    elseif bytesPerSamp == 2
                        nSamplesInd    = strcmp('num_EGF_samples',tempHeader(:,1));
                        nSamples       = sscanf(tempHeader{nSamplesInd,2},'%d');
                        %grab actual voltage data
                        tempData = fread(fid,nSamples,'int16'); %re-read as int16
                        obj.lfpData(1).lfpHighSamp{trialIterator}{j}  = (double(tempData)./2^15) .* obj.trialMetaData(trialIterator).lfp_scalemax(j); %voltages
                    end
                    fclose(fid);
                    
                    scanpix.fxchange.textprogressbar(j/length(lfp2load)*100);
                end
                if noEGFflag
                    scanpix.fxchange.textprogressbar('  NO DATA!');
                else
                    scanpix.fxchange.textprogressbar('  DONE!');
                end
            end
            
            obj.lfpData(1).lfpTet{trialIterator} = ceil(obj.trialMetaData(trialIterator).lfp_recordingChannel(:) ./ 4)';               
        end    
    end   
 %% This bit is from npix - makes waveforms
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
            %               pairs for params (see scanpix.npixUtils.getWaveforms)
            %
            % Outputs:
            %
            % See also: scanpix.npixUtils.getWaveforms
            %
            % LM 2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % deal with loading params
            if nargin == 1
                % use defaults, only add channel n from object metadata
                addParams = {'nch'; obj.trialMetaData(1).nChan};
            elseif strcmp(varargin{1},'ui')
                % UI dialoge
                prompts = {'N channels',                  'N waves / cluster', 'n chan / waveform', 'n samp / waveform', 'nSamplesPrePeak', 'apply CAR', 'unwhiten' };
                varargs = {'nch',                         'nwave',              'getnch',           'nsamp',             'prepeak',         'car',       'unwhite'  };
                defVals = { obj.trialMetaData(1).nChan,    250,                  5,                  40,                  0.375,             0,           0         };
                
                rtn = scanpix.helpers.makeCustomUIDialogue(prompts,defVals);
                if isempty(rtn)
                    warning('scaNpix::npix::loadWaves:Wave form loading aborted. That lacks class mate...');
                    return;
                end
                
%                 addParams = cell(2,size(rtn,1)-1);
%                 for i = 2:size(rtn,1)
%                     addParams{1,i-1} = varargs{i};
%                     addParams{2,i-1} = rtn{i,2};
%                 end
                addParams = cell(2,size(rtn,1));
                for i = 1:size(rtn,1)
                    addParams{1,i} = varargs{i};
                    addParams{2,i} = rtn{i,2};
                end
            else
                % name-value pairs
                addParams = {varargin{1:2:end};varargin{2:2:end}};
                % always add channel n unless already supplied
                if ~any(strcmp(addParams(1,:),'nch'))
                    addParams(:,end+1) = {'nch'; obj.trialMetaData(1).nChan};
                end

            end
            
            if any(strcmp(addParams(1,:),'save'))
                saveWFs = addParams{2,strcmp(addParams(1,:),'save')};
            else
                saveWFs = false;
            end

            for i = 1:length(obj.trialNames)
                [obj,tmpWF,tmpCH] = scanpix.npixUtils.extract_waveforms(obj,i,addParams{:});
                if saveWFs
                    waveforms = [tmpWF tmpCH];
                    save(fullfile(obj.dataPath{i},'waveforms.mat'),'waveforms');
                end
            end
        end
        
    end
  
end 