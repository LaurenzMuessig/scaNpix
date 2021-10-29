classdef dacq < handle
% DACQ data class. Let's you load in DACQ ephys data into a class
% object for analysis or data inspection
%
% Class constructor: scanpix.dacq
% Construct class obj  
%
% Usage:
%       obj = scanpix.dacq;
%       obj = scanpix.dacq(prmsMode);
%       obj = scanpix.dacq(prmsMode, uiFlag);
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
% LM 2020

    
    %% PROPERTIES
    properties % params %
        % Map object
        params                containers.Map
    end
    
    properties % meta data %
        dataPath(1,:)         string
        dataSetName           char
        trialNames(1,:)       string
        cell_ID(:,3)          double {mustBePositive, mustBeNonNan, mustBeNonzero}
    end
    
    properties % trial data %
        trialMetaData(1,:)    struct  % leave non-scalar as many fields and indexing into it, is not required very often
        posData               struct  = struct('XYraw',[],'XY',[],'direction',[],'speed',[],'linXY',[]);
        spikeData             struct  = struct('spk_Times',[],'spk_waveforms',[],'sampleT',[]);
        lfpData               struct  = struct('lfp',[],'lfpHighSamp',[],'lfpTet',[]);
    end
   
    properties % maps %
        maps                  struct  = struct('rate',[],'spike',[],'pos',[],'dir',[],'sACs',[],'OV',[],'speed',[]);
        linMaps               struct  = struct('linRate',[],'linPos',[],'linRateNormed',[]);
    end
    
%     properties(Dependent)
%         
%     end
    
    properties(Hidden)
        fileType              char    = '.set';
        uniqueCellIndex(:,1)  logical
        fields2spare          cell    = {'params','dataSetName','cell_ID'};  % spare this when deleting or rearranging data. make sure to add new properties that should be spared here!
        mapParams             struct  = scanpix.maps.defaultParamsRateMaps; 
        loadFlag              logical = false;                               % flag so we know something has been loaded into object
    end
    
    %% METHODS
    
    % constructor %
    methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = dacq(prmsMode,uiFlag)
            % dacq - creates class object 
            %
            % Syntax:  
            %       obj = dacq;
            %       obj = dacq(prmsMode);
            %       obj = dacq(prmsMode, uiFlag);
            %
            % Inputs:
            %    prmsMode - 'default' - uses default parameter (default)
            %               'ui'      - opens UI dialogue
            %               'file'    - load params from file
            %    uiFlag   - true (default)/false - skip UI set file
            %               selection dialogue if false
            %
            % Outputs:
            %    obj      - dacq object
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
    %%
    
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
            %    loadMode - cell array; {'all'} or any combination of {'pos','spikes','lfp'}
            %    varargin - list of set file names (ommit extensions)
            %
            % Outputs:
            %    
            % see also: obj.loadSet; obj.loadPos; obj.loadSpikes; obj.loadLFPs 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isempty(obj.dataPath)
                [obj.trialNames, obj.dataPath] = scanpix.helpers.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
                scanpix.helpers.selectTrials(obj.trialNames, obj);
            end
            
            if nargin < 2
                str = {'all','pos','spikes','lfp'};
                [select, loadCheck] = listdlg('PromptString','Select what data to load:','ListString',str,'ListSize',[160 100]);
                if ~loadCheck
                    warning('scaNpix: No data selected. Nothing is loaded. Boring...');
                    return;
                else
                    loadMode = str(select);
                end
            end
            
            if nargin < 3
                loadStr = obj.trialNames;
            else
                loadStr = obj.trialNames( ismember(obj.trialNames,varargin{1}) );   
            end
  
            cellTet = cell(length(loadStr),1);
            check0spikeCells = false; % flag
            trialInd = find(ismember(obj.trialNames,loadStr));
            
            % load trial data
            for i = 1:length(loadStr)
                
                obj.loadSet(trialInd(i)); % always load set file
                for j = 1:length(loadMode)
                    
                    switch lower(loadMode{j})
                        case 'all'
                            obj.loadPos(trialInd(i));
                            cellTet{i,1}    = obj.loadSpikes(trialInd(i));
                            obj.loadLFPs(trialInd(i));
                            check0spikeCells = true; % flag
                            break;
                        case 'pos'
                            obj.loadPos(trialInd(i));
                        case 'spikes'
                            cellTet{i,1}     = obj.loadSpikes(trialInd(i));
                            check0spikeCells = true; % flag
                        case 'lfp'
                            obj.loadLFPs(trialInd(i));
                        otherwise
                            ME = MException('scaNpix:load:invalidInput', ['' loadMode{j} ''' is not a valid data type selector. Try ''all'', ''pos'', ''spikes'' and/or ''lfp'' instead. ']);
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
                obj.cell_ID = [ allCellN ( 1:size(allCellN,1) )' ]; % also add a numerical index for cells (which is sometimes useful to have)
            end
            
%             scanpix.helpers.addMetaData(obj,obj.defaultMetaDataFields); % add the default meta data fields
            scanpix.helpers.preallocEmpty(obj,true,{'posData','spikeData','lfpData'});
            
            obj.loadFlag = true; % flag loading
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
            
            % legacy data isn't supported. multiplex pre-amp has different ADC range
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
            for i = 1:length(recordingChannel) %% There are occasional DACQs which only have EEG 1-4. This loop thereore needs to fail gracefully.
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
            setFile.trialType   = []; % no method yet for dacq data to gather that
            setFile.trackLength = []; % no method yet for dacq data to gather that
            
            % output
            if isempty(obj.trialMetaData)
                obj.trialMetaData = setFile;
            else
                obj.trialMetaData(trialIterator) = setFile;
            end
            
            if isempty(obj.dataSetName)
                obj.dataSetName = ['r_' obj.trialNames{trialIterator}(1:6)]; %%%%% THIS NEEDS FIX ONCE WE ARE CLEAR HOW TO HANDLE METADATA IN DACQ
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
                led_pos = floor( led_pos .* (obj.params('ScalePos2PPM') / str2double(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'))) );
                posHeader{strcmp(posHeader(:,1),'pixels_per_metre'),2} = num2str(obj.params('ScalePos2PPM'));
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
        
        %%
        function cell_ID = loadSpikes(obj, trialIterator)
            % loadSpikes - load tetrode data from DACQ files
            % Method for dacq class objects (hidden)
            % We will load spike times and waveforms for each cluster
            %
            % Syntax:  cell_ID = obj.loadSpikes(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %    cell_ID - numeric array of cell and tetrode IDs
            %
            % See also: dacq;
            %
            % TW/LM 2020 (adapted from org. SCAN function)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% To-Do
            % load MUA when no cut file is found
            %
            %%
            [~,nTetrodes] = scanpix.dacqUtils.findDACQFiles( [obj.dataPath{trialIterator} filesep], obj.trialNames{trialIterator}, 'tet' );
            if isempty(nTetrodes)
                ME = MException('scaNpix::loadSpikes:TetrodeFilesMissing', [' Can''t find any tetrode files for ' fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.set']) '.' ]);
                throw(ME);
            end
            
            clear scanpix.fxchange.textprogressbar  % get rid of persistent variable in 'textprogressbar' which can cause issues if you interrupt code during execution
            
            % initialise
            cell_ID                 = zeros(0, 2);
            [spikeTimes, waveforms] = deal({});
            scanpix.fxchange.textprogressbar(['Loading tetrode data for ' obj.trialNames{trialIterator} ' ']);
            
            % loop through recorded tetrodes
            for i = 1:length(nTetrodes)
                
                %%%%%%%%%%%%%% TETRODE FILE %%%%%%%%%%%%
                fid          = fopen(fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.' num2str(nTetrodes(i))]),'r','ieee-be');        % Open tet file. 'ieee-be' string is machine format, which needs to be 'big endian' to read times correctly. Default format doesnt work on my PC.
                data         = fread(fid,'int8'); % Read data (header and voltage samples). 'int8' converts each byte into a integer -126:126.
                
                % extract some info from header
                hdr          = char(abs(data(1:400)))';                % Translate header into text. Use 'abs' to avoid error message on negative integers after data starts.
                % n spikes
                nspk_ind     = strfind(hdr,'num_spikes ') + 11;        % Look for number of spikes marker
                nspk         = sscanf(hdr(nspk_ind:end), '%d');
                % Get sample rate (DSP timer Ticks/s)
                timebase_ind = strfind(hdr,'timebase ') + 9;
                timebase     = sscanf(hdr(timebase_ind:end), '%d');
                
                % read actual data
                ds           = strfind(hdr,'data_start') + 10;         % Look for data start marker
                % spike data
                spk_AmplRaw  = data( ds : ds + nspk*216 - 1);          % Take voltage sample values. (Time stamps still included)
                spk_AmplRaw  = reshape(spk_AmplRaw,54,4,nspk);         % Reshape into (sample,channel,spike)
                spk_AmplRaw  = spk_AmplRaw(5:end,:,:);                 % Cut off time stamp samples
                %%%%%%
                
                % time stamps
                fseek(fid,(ds-1),'bof');                               % Set file position indicator at data start (for reading time stamps, see below)
                spk_TimesRaw = fread(fid,nspk,'int32',212);            % read time stamps. 'nspk'=read this many samples, 'int32'=read 32bit integer, '212'=skip 212 bytes between reading values.
                spk_TimesRaw = spk_TimesRaw ./ timebase;
                %%%%%%
                fclose(fid);
                
                %%%%%%%%%%%%%% CUT FILE %%%%%%%%%%%%
                % read cut file
                if strcmp(obj.params('cutFileType'),'cut')
                    cutFileName = fullfile(obj.dataPath{trialIterator},[ obj.trialNames{trialIterator} obj.params('cutTag1') '_' num2str(nTetrodes(i)) obj.params('cutTag2') '.cut' ] );
                elseif strcmp(obj.params('cutFileType'),'clu')
                    cutFileName = fullfile(obj.dataPath{trialIterator}, [ obj.trialNames{trialIterator} obj.params('cutTag1') '.clu.' num2str(nTetrodes(i)) ] );
                else
                    ME = MException('scaNpix::loadSpikes:InvalidCutFileExtension', ['''' obj.params('cutFileType') '''' ' is not a valid extension for cut files - only ''cut'' or ''clu'' are allowed' ]);
                    throw(ME);
                end
                
                if ~exist(cutFileName,'file')
                    warning(['scaNpix: Can''t find ''' obj.trialNames{trialIterator} obj.params('cutTag1') '_' num2str(nTetrodes(i)) obj.params('cutTag2') '.cut''. No data loaded for Tet #' num2str(nTetrodes(i)) '.']);
                    continue;
                end
                
                fid = fopen(cutFileName,'r');
                
                if strcmp(obj.params('cutFileType'),'cut')
                    % get line with end of header (can be variable)
                    lineN = 0;
                    while 1
                        tline = fgetl(fid);                                 % Look for start of exact cut data.
                        if strncmp('Exact_cut_for',tline,12), lineN = lineN+1; break, end      %
                        lineN = lineN+1;                                              % The final returned tline is the header for this.
                    end
                else
                    lineN = 1; % n.b. first line in clu files is n of clusters found
                end
                %
                frewind(fid);
                cut = cell2mat( textscan(fid,'%n','headerlines',lineN,'delimiter','/b') );
                fclose(fid);
                
                % sanity check
                if length(spk_TimesRaw) ~= length(cut)
                    warning(['scaNpix: Trial ', obj.trialNames{trialIterator}, ': Number of spikes in cut and tetrode file for tetrode ' num2str(nTetrodes(i)) ' don''t match. Cut file corrupt?. Data loading for this tetrode will be skipped']);
                    continue;
                end
                
                % remove spikes in overhang
                overhangInd  = spk_TimesRaw < obj.trialMetaData(trialIterator).duration;
                cut          = cut(overhangInd); % remove spikes in overhang
                spk_TimesRaw = spk_TimesRaw(overhangInd);
                spk_AmplRaw  = spk_AmplRaw(:,:,overhangInd);
                
                clustIds = unique(cut(cut~=0));  %  Note that .clu files have no convention of a cluster 0 that contains unsorted spikes (i.e. MUA and noise) - maybe we should make one?
                
                if isempty(clustIds)
                    scanpix.fxchange.textprogressbar(i/length(nTetrodes)*100);
                    continue;
                end
                
                [cut, sortInd]           = sort(cut);               % this will ensure output of accumarray is sorted properly (spike times in each cluster)
                spk_TimesRaw             = spk_TimesRaw(sortInd);   % sort other arrays as well accordingly
                spk_AmplRaw              = spk_AmplRaw(:,:,sortInd);
                % get spike times for each cluster
                temp_ST                  = accumarray(cut(cut~=0), spk_TimesRaw(cut~=0,1), [max(unique(cut)) 1], @(x) {x});
                ind_empty                = ~cellfun('isempty',temp_ST);
                temp_ST                  = temp_ST(ind_empty);
                % bit more complicated for waveform amplitudes
                % first convert to uV...
                chGains                  = repmat(obj.trialMetaData(trialIterator).fullscale(nTetrodes(i),:),[50,1,sum(overhangInd)]);
                spk_AmplRaw1             = spk_AmplRaw ./128 .* chGains;
                % ... then assign to clusters
                temp                     = shiftdim(repmat(cut, 1, 50, 4), 1); % all cluster IDs for raw amplitude array
                temp                     = temp(:);
                temp2                    = spk_AmplRaw1(:);
                temp_SA                  = accumarray(temp(temp~=0), temp2(temp~=0), [max(unique(temp)) 1],@(x) {x});
                temp_SA                  = cellfun(@(x) shiftdim( reshape(x,50,4,[] ), 2),temp_SA, 'UniformOutput',0); % restore old format ('reshape') & bring into same format as 'temp_ST' ('shiftdim')
                temp_SA                  = temp_SA(ind_empty);
                % accumulate output data
                spikeTimes(end+1:end+length(temp_ST),1) = temp_ST;
                waveforms(end+1:end+length(temp_ST),1)  = temp_SA;
                cell_ID(end+1:end+length(temp_ST),:)    = [clustIds, nTetrodes(i) * ones(length(clustIds),1)]; % [cell n, tet n]
                
                scanpix.fxchange.textprogressbar(i/length(nTetrodes)*100);
            end
            obj.spikeData(1).spk_Times{trialIterator}     = spikeTimes;
            obj.spikeData(1).spk_waveforms{trialIterator} = waveforms;
            obj.spikeData(1).sampleT{trialIterator}       = [];  % not relevant for dacq data, but we need to keep format consistent with neuropix data  
            
            scanpix.fxchange.textprogressbar('  DONE!');           
        end
        
        %%
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
end
