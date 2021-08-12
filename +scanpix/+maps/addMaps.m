function addMaps(obj, mapType, trialInd, varargin )
% addMaps - create different types of firing rate maps for data in object
% package: scanpix.maps
%
% We can add spatial rate maps ('rate'), directional rate maps ('dir') or
% linearised rate maps ('lin') in case of linear track data. This function 
% is mainly to allow for flexibly choosing which parameters should be used 
% for the map construction. 
% For the parameter space check 'scanpix.maps.defaultParamsRateMaps'. 
%
% Syntax:
%       scanpix.maps.addMaps(obj)
%       scanpix.maps.addMaps(obj,mapType)
%       scanpix.maps.addMaps(obj,mapType, trialInd)
%       scanpix.maps.addMaps(obj,mapType, trialInd, varargin )
%       scanpix.maps.addMaps(obj,mapType, [], varargin )
%       scanpix.maps.addMaps(__, 'load' )
%       scanpix.maps.addMaps(__, 'ui' )
%       scanpix.maps.addMaps(__, prmsStruct )
%       scanpix.maps.addMaps(__, Name-Value comma separated list )
%
% Inputs:
%    obj      - dacq or npix class object
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
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% quite a lot of input parsing/checking to do here %%%%%%%%
if ~obj.loadFlag
    warning('scaNpix::maps::addMaps: You need to load some data first before you can make any maps... Fairly obvious if you ask me.');
    return;
end

if nargin < 2
    str = {'rate','dir','lin','objVect','speed'};
    [select, loadCheck] = listdlg('PromptString','Select what maps to make:','ListString',str,'ListSize',[160 100],'SelectionMode','Single');
    if ~loadCheck
        warning('scaNpix::maps::addMaps: No data selected. No maps will be created. Boring...');
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
        prms = obj.mapParams;
    else
        % use user defined params
        classFolder = what('scanpix');
        try
            temp    = load( fullfile(classFolder.path,'files', obj.params('myRateMapParams') ) ); %
        catch
            ME = MException('scaNpix::maps::addMaps:''rateMapParamsNotFound', ['' fullfile(classFolder.path,'files', obj.params('myRateMapParams') ) ''' doesn''t exist buddy... ']);
            throw(ME);
        end
        f            = fieldnames(temp);
        prms         = temp.(f{1});
        warning(['scaNpix::maps::addMaps: Using user-defined parameters to make maps loaded from ' fullfile(classFolder.path,'files', obj.params('myRateMapParams') ) '.'] );
    end
end

if nargin == 4 && ischar(varargin{1})
    if strcmpi(varargin{1},'default')
        % use defaults
        prms = scanpix.maps.defaultParamsRateMaps;
        
    elseif strcmpi(varargin{1},'load')
        % load from file
        [fName, dataDir] = uigetfile(fullfile(classFolder.path,'files', '*.mat'), 'Select map params to load.');
        % fail gracefully
        if isnumeric(fName)
            warning('scaNpix:maps:addMaps: Loading Map Params aborted. Will use defaults instead!');
            prms         = obj.mapParams;
        else
            temp         = load( fullfile(dataDir, fName) );
            f            = fieldnames(temp);
            prms         = temp.(f{1});
        end
        
    elseif strcmpi(varargin{1},'ui')
        % open UI dialogue to grab parameters
        prms             = obj.mapParams;
        prompts          = fieldnames(prms);
        defaultVals      = struct2cell(prms);
        output           = scanpix.helpers.makeCustomUIDialogue(prompts, defaultVals);
        % exit gracefully
        if isempty(output)
            warning('scaNpix::maps::addMaps: Aborted changing parameters for rate map generation - will use those currently set in object.');
            prms         = obj.mapParams;
        else
            % format as structure
            prms         = cell2struct(output(:,2), output(:,1), 1);
        end
        
    end
end

% Name-Value pair or struct input to change some params
% explicitly
if ~isempty(varargin) && (nargin > 4 || isstruct(varargin{1}) )
    prms = obj.mapParams;
    if ischar(varargin{1})                                                           %
        for i=1:2:length(varargin);   prms.(varargin{i}) = varargin{i+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for i=1:length(f);   prms.(f{i}) = s.(f{i});   end                        %
    end
end

% these params get set from general data params container
prms.PosFs               = obj.params('posFs');
% prms.ppm                 = obj.params('ppm');
prms.speedFilterLimits   = [prms.speedFilterLimitLow prms.speedFilterLimitHigh];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make some maps
switch lower(mapType)
    case 'rate'
        for i = trialInd
            [ obj.maps(1).rate{i}, obj.maps(1).pos{i}, obj.maps(1).spike{i} ] = scanpix.maps.makeRateMaps(obj.spikeData.spk_Times{i}, obj.posData.XY{i}, obj.spikeData.sampleT{i}, obj.trialMetaData(i).ppm, obj.posData.speed{i}, prms );
        end
        
        % pad maps so size is the same for all maps is set
        mapSz = reshape(cell2mat(cellfun(@(x) size(x{1}),obj.maps.rate,'uni',0)),2,[]);
        maxSz = max(mapSz,[],2);
        padSz = ceil((maxSz-mapSz) ./ 2);
        
        mapType = {'rate','pos','spike'};
        for i = 1:length(mapType)
            for j = trialInd
                obj.maps.(mapType{i}){j} = cellfun(@(x) padarray(x,[padSz(1,j) padSz(2,j)],NaN), obj.maps.(mapType{i}){j},'uni',0);
                if any(mapSz(:,j) + 2*padSz(:,j) > maxSz)
                    obj.maps.(mapType{i}){j} = cellfun(@(x) x(1:maxSz(1),1:maxSz(2)),obj.maps.(mapType{i}){j},'uni',0); % might need to trim off one row/column
                end
            end
        end
        
    case 'dir'
        for i = trialInd
            obj.maps(1).dir{i} = scanpix.maps.makeDirMaps( obj.spikeData.spk_Times{i}, obj.posData.direction{i}, obj.spikeData.sampleT{i}, obj.posData.speed{i},  prms  );
        end
        
    case 'lin'
        % map making a bit more involved mostly because we need some extra parameters we can't infer from data alone - should check if we
        % can solve this more elegantly
        skipNextUI = false;
        for i = trialInd
            % we need the type of the track (for how smoothing is done)
            if ~isempty(obj.trialMetaData(trialInd).trialType)
                trackProps.type = obj.trialMetaData(trialInd).trialType;
            else
                uiInput = inputdlg({'linear track type', 'linear track length (cm)'},'',1,{'sqtrack','62.5'});
                if isempty(uiInput)
                    warning('scaNpix::maps::addMaps:No track properties, no linear rate maps...');
                    return;
                else
                    trackProps.type = uiInput{1};
                    obj.trialMetaData(trialInd).trialType = uiInput{1};
                    trackProps.length = str2double(uiInput{2});
                    obj.trialMetaData(trialInd).trackLength = str2double(uiInput{2});
                    skipNextUI = true;
                end
            end
            
            % we need the length of the track for the pos scaling to work
            if ~isempty(obj.trialMetaData(trialInd).trackLength)
                trackProps.length= obj.trialMetaData(trialInd).trackLength;
            elseif ~skipNextUI
                uiInput = inputdlg({'linear track length (cm)'},'',1,{'62.5'});
                if isempty(uiInput)
                    warning('scaNpix::maps::addMaps:No track length, no linear rate maps...');
                    return;
                else
                    trackProps.length = str2double(uiInput{1});
                    obj.trialMetaData(trialInd).trackLength = str2double(uiInput{1});
                end
            end
            trackProps.ppm   = obj.trialMetaData(i).ppm;
            trackProps.posFs = obj.params('posFs'); 
            [obj.linMaps(1).linRate{i},obj.linMaps(1).linPos{i},obj.posData(1).linXY{i},obj.linMaps(1).linRateNormed{i}] = scanpix.maps.makeLinRMaps(obj.spikeData.spk_Times{i}, obj.posData.XY{i},...
                                                                                                                                obj.spikeData.sampleT{i}, obj.posData.direction{i},obj.posData.speed{i}, trackProps, prms );
        end
        
    case 'objvect'
        for i = trialInd
            if ~isempty(obj.trialMetaData(i).objectPos)
                obj.maps(1).OV{i} = scanpix.maps.makeOVMap( obj.spikeData.spk_Times{i}, obj.posData.XY{i}, obj.spikeData.sampleT{i}, obj.trialMetaData(i).objectPos, obj.trialMetaData(i).ppm,  prms  );
            end
        end
        
    case 'speed'
        
        for i = trialInd
            obj.maps(1).speed{i} = scanpix.maps.makeSpeedMap( obj.spikeData.spk_Times{i}, obj.posData.speed{i}, obj.trialMetaData(i).duration,  prms  );
        end
        
    otherwise
        ME = MException('scaNpix::maps::addMaps:invalidMapType', ['' mapType ''' is not yet a supported map type. You need to invent that one yourself I am afraid... ']);
        throw(ME);
end

end