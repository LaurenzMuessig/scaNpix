classdef npix < handle
    % npix data class. Let's you load in neuropixel data into a class
    % object for analysis or data inspection
    %
    % Class constructor: scanpix.npix
    % Construct class obj
    %
    % Usage:
    %       obj = scanpix.npix;
    %       obj = scanpix.npix(prmsMode);
    %       obj = scanpix.npix(prmsMode, uiFlag);
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
        chanMap               struct
    end
    
    properties % meta data %
        dataPath(1,:)         string
        dataSetName(1,:)      char
        trialNames(1,:)       string
        cell_ID(:,3)          double %{mustBePositive, mustBeNonNan, mustBeNonzero}
        cell_Label(:,1)       string
    end
    
    properties % trial data %
        trialMetaData(1,:)    struct
        posData               struct  = struct('XYraw',[],'XY',[],'direction',[],'speed',[],'linXY',[],'sampleT',[]);
        spikeData             struct  = struct('spk_Times',[],'spk_waveforms',[],'sampleT',[]);  %%% DEAL WITH WAVEFORMS?
        lfpData               struct  = struct('lfp',[]);
    end
    
    properties % maps %
        maps                  struct  = struct('rate',[],'spike',[],'pos',[],'dir',[],'sACs',[],'OV',[],'speed',[]);
        linMaps               struct  = struct('linRate',[],'linPos',[],'linRateNormed',[]);
    end
    
    %     properties(Dependent)
    %
    %     end
    
    properties(Hidden)
        fileType              char    = '.ap.bin';
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
        function obj = npix(prmsMode,uiFlag)
            % dacq - creates class object
            %
            % Syntax:
            %       obj = npix;
            %       obj = npix(prmsMode);
            %       obj = npix(prmsMode, uiFlag);
            %
            % Inputs:
            %    prmsMode - 'default' - uses default parameter (default)
            %               'ui'      - opens UI dialogue
            %               'file'    - load params from file
            %    uiFlag   - true (default)/false - skip UI set file
            %               selection dialogue if false
            %
            % Outputs:
            %    obj      - npix object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            if nargin == 0
                obj.params = scanpix.helpers.getParams(obj,'default');
                uiFlag     = true;
            else
                obj.params = scanpix.helpers.getParams(obj,prmsMode);
                if nargin < 2
                    uiFlag = true;
                end
            end
            
            if isempty(obj.params)
                warning('scanNpix: No params selected. Using defaults.');
                obj.params = scanpix.helpers.defaultParamsContainer(class(obj));
            end
            
            if uiFlag
                [obj.trialNames, obj.dataPath] = scanpix.helpers.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
                if isempty(obj.trialNames); return; end
                scanpix.helpers.selectTrials(obj.trialNames, obj);
            end
            
            if ~isempty(obj.dataPath)
                obj.params('defaultDir') =[fileparts(obj.dataPath{1}) filesep];  % use parent directory rather than absolut as this makes more sense for multiple trial data
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
                
                obj.loadMeta(trialInd(i)); % load some meta data
                               
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
    methods(Hidden)
        %%
        function loadMeta(obj, trialIterator)
            % loadMeta - load meta data for neuropixel data
            % Method for scanpix class objects (hidden)
            %
            % Syntax:  obj.loadMeta(trialIterator)
            %
            % Inputs:
            %    trialIterator - numeric index for trial to be loaded
            %
            % Outputs:
            %
            % See also: scanpix;
            %
            % LM 2021
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % first load xml file
            metaXMLFileInfo = dir(fullfile(obj.dataPath{trialIterator},'*.xml'));
            metaXMLFile     = scanpix.fxchange.xml_read(fullfile(obj.dataPath{trialIterator},metaXMLFileInfo.name)); %
            
            f = fieldnames(metaXMLFile);
            for i = 1:length(f)
                if ~strcmpi(f{i},'comment')
                    obj.trialMetaData(trialIterator).(f{i}) = metaXMLFile.(f{i});
                end
            end
            % reshape into more convenient format ([minX maxX; minY maxY])
            obj.trialMetaData(trialIterator).envBorderCoords = [min(obj.trialMetaData(trialIterator).envBorderCoords([1,3])),max(obj.trialMetaData(trialIterator).envBorderCoords([1,3])); min(obj.trialMetaData(trialIterator).envBorderCoords([2,4])),max(obj.trialMetaData(trialIterator).envBorderCoords([2,4]))];

            % legacy: older xml files won't contain 'objectPos' field
            if ~isfield(metaXMLFile,'objectPos')
                obj.trialMetaData(trialIterator).objectPos = [];
            end
            %
            if ~isfield(metaXMLFile,'posFs')
                obj.trialMetaData(trialIterator).posFs = 50; % HARCODED ATM! Should maybe be added to xml file?
            end
            
            obj.trialMetaData(trialIterator).ppm = [];
            obj.trialMetaData(trialIterator).ppm_org = [];
            obj.trialMetaData(trialIterator).trackLength = []; % add field to xml?
            % spikeGLX meta data %
            metaDataFile = dir(fullfile(obj.dataPath{trialIterator},'*.ap.meta'));
            fidMeta  = fopen(fullfile(metaDataFile.folder,metaDataFile.name),'r');
            C        = textscan(fidMeta, '%[^=] = %[^\r\n]');
            fclose(fidMeta);
%             obj.trialMetaData(trialIterator).nChanTot = sscanf(C{2}{strcmp(C{1},'nSavedChans')},'%d');
%             obj.trialMetaData(trialIterator).nChanAP  = sscanf(C{2}{strcmp(C{1},'snsApLfSy')},'%d%*%*');
            %%%% Do we want to add more info from metafile?? %%%%%%%%%%%%%%%%%
            
            % load channel map
            chanMapInfo = dir(fullfile(obj.dataPath{trialIterator},'*ChanMap*'));
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
            
            if isempty(obj.dataSetName)
                obj.dataSetName = ['r' num2str(metaXMLFile.animal) '_' num2str(metaXMLFile.date)];
            end
            

        end
        
        %%
        function loadPos(obj, trialIterator)
            % loadPos - load position data for neuropixel data
            % Method for scanpix class objects (hidden)
            %
            % Syntax:  obj.loadPos(trialIterator)
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
            
            fprintf('Loading pos data for %s .......... ', obj.trialNames{trialIterator});
            
            
            %% process data
            
            % open pos data the format is [frame count, greenXY, redXY winSzX, winSzY, timeStamp possibly other Data ]
            fName = dir(fullfile(obj.dataPath{trialIterator},'trackingData', '*.csv'));
            if isempty(fName)
                disp(['Can''t find csv file in ' obj.dataPath{trialIterator} '. Come on mate.']);
                return;
            end
            
            fID = fopen(fullfile(fName.folder,fName.name),'rt');
            csvData = textscan(fID,'%u%f%f%f%f%u%u%f','HeaderLines',1,'delimiter',',');
            fclose(fID);
            
            % led data - same format as for dacq
            if strcmp(obj.trialMetaData(trialIterator).LEDfront,'green')
                led          = [csvData{2}, csvData{3}]; % xy coords
                led(:,:,2)   = [csvData{4}, csvData{5}]; % xy coords
            else
                led          = [csvData{4}, csvData{5}]; % xy coords
                led(:,:,2)   = [csvData{2}, csvData{3}]; % xy coords
            end
            led(led==0)      = NaN;
            
            % sample Times
            timeStamps       = csvData{8};
            sampleT          = scanpix.npixUtils.convertPointGreyCamTimeStamps(timeStamps); % starts @ 0
            
            % in case logging point grey data was corrupt
            if all(sampleT == 0)
                sampleT          = (0:1/obj.trialMetaData(trialIterator).posFs:length(led)/obj.trialMetaData(trialIterator).posFs)';
                sampleT          = sampleT(1:length(led)); % pretend we have perfect sampling
                frameCount       = 1:length(led); % pretend we are not missing any frames
                obj.trialMetaData(trialIterator).BonsaiCorruptFlag = true;
                warning('Point Grey data corrupt!');
            else
                frameCount       = csvData{1} - csvData{1}(1) + 1;
                obj.trialMetaData(trialIterator).BonsaiCorruptFlag = false;
            end
            
            % deal with missing frames (if any) - this currently doesn't take into
            % account if 1st frame(s) would be missing, but I am not sure this would
            % actually ever happen (as 1st fame should always be triggered fine)
            % first check if there are any...
            missFrames       = find(~ismember(1:length(frameCount),frameCount));
            nMissFrames      = length(missFrames);
            if ~isempty(missFrames)
                fprintf('Note: There are %i missing frames in tracking data for %s.\n', nMissFrames, obj.trialMetaData(trialIterator).filename);
                
                temp                   = zeros(length(led)+nMissFrames, 2, obj.trialMetaData(trialIterator).nLEDs);
                temp(missFrames,:,:)   = nan;
                temp(temp(:,1)==0,:,:) = led;
                led                    = temp;
                
                % interpolate sample times
                interp_sampleT         = interp1(double(frameCount), sampleT, missFrames);
                temp                   = zeros(length(led),1);
                temp(missFrames,1)     = interp_sampleT;
                temp(temp(:,1) == 0,1) = sampleT;
                sampleT                = temp;
            end
            
            ppm = nan(2,1);
            % estimate ppm
            if isempty(obj.trialMetaData(trialIterator).envBorderCoords)
                winSz  = [double(csvData{6}(1)) double(csvData{7}(1))];
                ppm(:) = mean(winSz ./ (obj.trialMetaData(trialIterator).envSize ./ 100) );
            else
                % this case should be default
                winSz  = [abs(obj.trialMetaData(trialIterator).envBorderCoords(1)-obj.trialMetaData(trialIterator).envBorderCoords(3)), abs(obj.trialMetaData(trialIterator).envBorderCoords(2)-obj.trialMetaData(trialIterator).envBorderCoords(4))];
                ppm(:) = round( mean( winSz ./ (obj.trialMetaData(trialIterator).envSize ./ 100) ) );
            end
            
            %% post process - basically as scanpix.dacqUtils.postprocess_data_v2
            % scale data to standard ppm if desired
            if ~isempty(obj.params('ScalePos2PPM'))
                scaleFact = (obj.params('ScalePos2PPM')/ppm(1));
                led = floor(led .* scaleFact);
                ppm(1) = obj.params('ScalePos2PPM');
                obj.trialMetaData(trialIterator).objectPos = obj.trialMetaData(trialIterator).objectPos .* scaleFact;
                obj.trialMetaData(trialIterator).envBorderCoords = obj.trialMetaData(trialIterator).envBorderCoords .* scaleFact;
            end
            
            % remove tracking errors (i.e. too fast)
            for i = 1:2
                % speed
                pathDists        = sqrt( diff(led(:,1,i),[],1).^2 + diff(led(:,2,i),[],1).^2 ) ./ ppm(1); % % distances in m
                tempSpeed        = pathDists ./ diff(sampleT); % cm/s
                tempSpeed(end+1) = tempSpeed(end);
                speedInd = tempSpeed > obj.params('posMaxSpeed');
                % env borders
                if isempty(regexp(obj.trialMetaData(trialIterator).trialType,'circle','once'))
                	envSzInd = led(:,1,i) < 0.95*obj.trialMetaData(trialIterator).envBorderCoords(1,1) | led(:,1,i) > 1.05*obj.trialMetaData(trialIterator).envBorderCoords(1,2) | led(:,2,i) < 0.95*obj.trialMetaData(trialIterator).envBorderCoords(2,1) | led(:,2,i) > 1.05*obj.trialMetaData(trialIterator).envBorderCoords(2,2); 
                else
                    %%%% NEED TO ADD CIRCLE FIT %%%%
                end
                % filter out
                led(speedInd | envSzInd,:,i) = NaN;
            end
            
            % interpolate missing positions
            for i = 1:2
                missing_pos = find(isnan(led(:,1,i)));
                ok_pos      = find(~isnan(led(:,1,i)));
                for j = 1:2
                    led(missing_pos, j, i)                            = interp1(ok_pos, led(ok_pos, j, i), missing_pos, 'linear');
                    led(missing_pos(missing_pos > max(ok_pos)), j, i) = led( max(ok_pos), j, i);
                    led(missing_pos(missing_pos < min(ok_pos)), j, i) = led( min(ok_pos), j, i);
                end
            end
            
            % smooth 
            kernel         = ones( ceil(obj.params('posSmooth') * obj.params('posFs')), 1)./ ceil( obj.params('posSmooth') * obj.params('posFs') ); % as per Ephys standard - 400ms boxcar filter
            % Smooth lights individually, then get direction.
            smLightFront   = imfilter(led(:, :, 1), kernel, 'replicate');
            smLightBack    = imfilter(led(:, :, 2), kernel, 'replicate');
            
            correction                              = obj.trialMetaData(trialIterator).LEDorientation(1); %To correct for light pos relative to rat subtract angle of large light
            dirData                                 = mod((180/pi) * ( atan2(smLightFront(:,2)-smLightBack(:,2), smLightFront(:,1)-smLightBack(:,1)) ) - correction, 360); % 
            obj.posData(1).direction{trialIterator} = dirData(:);
            % Get position from smoothed individual lights %%
            wghtLightFront = 1-obj.params('posHead');
            wghtLightBack  = obj.params('posHead');
            xy = (smLightFront .* wghtLightFront + smLightBack .* wghtLightBack);  %
            
            % pos data
            obj.posData(1).XYraw{trialIterator}        = xy;
            obj.posData(1).XY{trialIterator}           = [double( floor(xy(:,1)) + 1 ), double( floor(xy(:,2)) + 1 )];
            obj.posData(1).sampleT{trialIterator}      = sampleT; % this is redundant (I think)
            
            obj.trialMetaData(trialIterator).ppm       = ppm(1);
            obj.trialMetaData(trialIterator).ppm_org   = ppm(2);
            
            % scale position
            boxExt = obj.trialMetaData(trialIterator).envSize / 100 * obj.trialMetaData(trialIterator).ppm;
            scanpix.maps.scalePosition(obj, trialIterator, 'envDimensions',boxExt);
            
            % running speed
            pathDists                                  = sqrt( (obj.posData(1).XY{trialIterator}(1:end-1,1) - obj.posData(1).XY{trialIterator}(2:end,1)).^2 + (obj.posData(1).XY{trialIterator}(1:end-1,2) - obj.posData(1).XY{trialIterator}(2:end,2)).^2 ) ./ ppm(1) .* 100; % distances in cm
            obj.posData(1).speed{trialIterator}        = pathDists ./ diff(sampleT); % cm/s
            obj.posData(1).speed{trialIterator}(end+1) = obj.posData(1).speed{trialIterator}(end);

            fprintf('  DONE!\n');
            
        end
        
        %%
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
            % load spike times
            spikeTimes         = readNPY(fullfile(obj.dataPath{trialIterator},'spike_times.npy'));
            spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
            obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
            % load cluster IDs
            clustIDs           = readNPY(fullfile(obj.dataPath{trialIterator},'spike_clusters.npy')) + 1; % 0 based index
            % now we need to remove bad clusters and spike times outside trial
            unClustIDs         = unique(clustIDs);
            if loadFromPhy    
                % phy curated output - only load clusters labelled good
                clu_info       = tdfread(fullfile(obj.dataPath{trialIterator},'cluster_info.tsv'),'tab');
                goodLabel      = strcmp(cellstr(clu_info.group),'good');
                good_clusts    = clu_info.id(goodLabel) + 1;
                clu_Depth      = clu_info.depth(goodLabel);  % I think this amd the following seem to not give correct results in somce cases?? Phy Issue?? 
                clu_Ch         = clu_info.ch(goodLabel) + 1; % this is 0 based 
                cluLabel       = string(clu_info.group);
                cluLabel       = cluLabel(goodLabel);
            else
                % raw kilosort output - we'll take everything in this case and can
                % filter later if needed
                ks_labels      = tdfread(fullfile(obj.dataPath{trialIterator},'cluster_KSLabel.tsv'),'tab'); % we need the raw kilosort output here before you ran Phy as this file gets overwritten! - maybe point to backup folder?
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
            good_clusts          = good_clusts(indSort);
            cluLabel             = cluLabel(indSort);
            clu_Ch               = clu_Ch(indSort);
            
            
            %% output
            obj.spikeData(1).spk_Times{trialIterator} = spikeTimesFin;
            endIdxNPix                                = min( [ length(syncTTLs), find(syncTTLs < obj.trialMetaData(trialIterator).duration + syncTTLs(1),1,'last') + 1]);
            obj.spikeData(1).sampleT{trialIterator}   = syncTTLs(1:endIdxNPix) - syncTTLs(1);
            % we'll only need to grab this once - this assumes data was clustered together, but wouldn't make much sense otherwise?
            % what about reload?
            if trialIterator == 1
                obj.cell_ID    = [good_clusts, clu_Depth, clu_Ch];
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
