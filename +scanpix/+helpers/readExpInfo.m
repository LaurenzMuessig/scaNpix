function [ expInfo, maxNTrial ] = readExpInfo( cribSheetPath, varargin )
% readExpInfo - Read in trial meta data from a crib sheet to be able to e.g. load data in big analysis loop or load specific data into Matlab. 
% We assume the following:
%   - the sheet is an excel file and it needs to have the following format 
%     (column headers): 
%     - filePath (full path to data) 
%     - animal(animal ID as 'rXXX' or 'mXXX') 
%     - trialName ('setFileName' w/o extension)
%     - experiment_ID (numeric index that uniquely identifies experiments within animals)
% 
% Note that when no arguments are passed we will open a UI dialogue to find sheet on disk
%
% Usage:    expInfo = readCribsheet
%           expInfo = readCribsheet( cribSheetPath );
%           expInfo = readCribsheet( cribSheetPath, method )
%           expInfo = readCribsheet( cribSheetPath, method, 'inputName', inputVal, .. etc .. )
%
% Inputs:   cribSheetPath - full path to file (including name & extension) 
%           method        - string; defines how data is read; 
%                           'single' (read each line of sheet individually, i.e. no grouping) - default
%                           'singleTrial' (read 1 specific trial - needs additional params to work)
%                           'exp' (read full sheet and group by animal + experiment_ID) 
%                           'singleExp' (read 1 specific experiment - needs additional params to work)
%           varargin      - name-value: comma separated list of name-value pairs
%
%
% Outputs:  
%           expInfo   - structure with meta data 
%           maxNTrial - max trial number in any experiment (good to know for pre-allocation) 
%
%
% LM 2020 / 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% machine/format specific fields %%%
defaultSheetPath   = 'S:\1postDoc\Neuropixels\rawData\masterSheetNpixData.xlsx';
%%%

standardFieldNames = {'filePath','animal','trialName','experiment_ID'}; % these fields are standard and need to be included on sheet (except for 'notes')

%% params
defaultMethod       = 'single';
defaultSheet        = 1; 
defaultAnimalID     = ''; 
defaultExperimentID = [];
defaultTrialName    = '';

p = inputParser;
validMethodInput = @(x) any(strcmp(x,{'single','singletrial','exp','singleexp'}));
addOptional(p,'method', defaultMethod, validMethodInput);
addParameter(p,'sheetN', defaultSheet);
addParameter(p,'animal', defaultAnimalID, @ischar);
addParameter(p,'experiment', defaultExperimentID);
addParameter(p,'trialname', defaultTrialName, @ischar);

parse(p,varargin{:});


%% parse input
if nargin == 0
    uiAnswer = inputdlg({'Path to Cribsheet:', 'load method', 'sheetN', 'experiment_ID', 'animal'},'Gather data for Data Loading',[1 100],{defaultSheetPath,defaultMethod,num2str(defaultSheet),'',''});
    
    if isempty(uiAnswer)
        warning('Loading. This code terminates here.');
        expInfo   = [];
        maxNTrial = [];
        return
    end
    
    cribSheetPath = uiAnswer{1};
    method        = uiAnswer{2};
    sheetN        = str2double(uiAnswer{3});
    expID         = str2double(uiAnswer{4});
    animal        = str2double(uiAnswer{5});
else
    method        = p.Results.method;
    sheetN        = p.Results.sheetN;
    expID         = p.Results.experiment;
    animal        = p.Results.animal;
end

% sanity check
if exist(cribSheetPath,'file') ~= 2
    error(['Can''t find crib sheet for metadata in ''' cribSheetPath '''. You have to do better than that.']);
end

%% read in data
cribSheet = readtable( cribSheetPath,'Sheet', sheetN ); % read cribsheet
% you might leave trial name empty in case data folder only contains data
% for a single trial
if ~iscell(cribSheet.trialName)
    cribSheet.trialName = cell(height(cribSheet),1);
end

% get a list of other fields in cribsheet - add these to output as well
addFieldNames = cribSheet.Properties.VariableNames( ~ismember(lower(cribSheet.Properties.VariableNames), lower(standardFieldNames)) );

anNum = sscanf([cribSheet.animal{:}],'%*c%d');  % a numerical animal identifier is also helpful - need to remove the 'r' or 'm'. 

switch lower(method)
    
    case 'single'
        % read every trial in sheet seperately
        expInfo.animal       = cribSheet.animal;
        expInfo.anNum        = anNum;
        expInfo.Tnames       = cribSheet.trialName;
        expInfo.fullPath     = fullfile(cribSheet.filePath,cribSheet.trialName);
        % additional data (meta data)
        if ~isempty(addFieldNames)
            for j = 1:length(addFieldNames)
                if ~iscell(cribSheet.(addFieldNames{j}))
                    expInfo.(addFieldNames{j}) = num2cell(cribSheet.(addFieldNames{j}), 2);
                else
                    expInfo.(addFieldNames{j}) = cellfun(@(x){x}, cribSheet.(addFieldNames{j}),'uni',0);
                end
            end
        end
        
        maxNTrial = 1;
        
    case 'singletrial'
        % read in single trial
        ind4trial            = strcmp(cribSheet.animal,animal) & strcmp(cribSheet.trialName,p.Results.trialname);
        expInfo.animal       = cribSheet.animal(ind4trial);
        expInfo.anNum        = anNum(ind4trial);
        expInfo.Tnames       = cribSheet.trialName(ind4trial);
        expInfo.fullPath     = fullfile(cribSheet.filePath(ind4trial),cribSheet.trialName(ind4trial));
        % additional data (meta data)
        if ~isempty(addFieldNames)
            for j = 1:length(addFieldNames)
                expInfo.(addFieldNames{j}){1} = cribSheet.(addFieldNames{j})(ind4trial);
            end
        end
        
        maxNTrial = 1;

    case 'exp'  
        % group data by experiment
        [~,~,expUIDInd] = unique( [anNum, cribSheet.experiment_ID], 'rows' );  % 'expUIDInd' is a unique identifier for experiments (as the field 'experiment_ID' in the sheet is only unique within animal).
        uniqueExpList   = unique(expUIDInd);
        
        for i = 1:length(uniqueExpList)   % i=iterator for EXPERIMENT
            trialIndForExp             = find( expUIDInd == uniqueExpList(i) );
            expInfo.animal{i,1}        = cribSheet.animal{trialIndForExp(1)};
            expInfo.anNum{i,1}         = anNum(trialIndForExp);
            expInfo.Tnames{i,1}        = cribSheet.trialName( trialIndForExp );  
            expInfo.fullPath{i,1}      = fullfile(cribSheet.filePath(trialIndForExp),cribSheet.trialName(trialIndForExp)); %%%%%%%
            expInfo.experiment_ID{i,1} = cribSheet.experiment_ID( trialIndForExp ); 
            
            % additional data (meta data)
            if ~isempty(addFieldNames)
                for j = 1:length(addFieldNames)
                    expInfo.(addFieldNames{j}){i,1} = cribSheet.(addFieldNames{j})(trialIndForExp);
                end    
            end 
            
        end
        maxNTrial = max( accumarray(expUIDInd,1) );
        
    case 'singleexp'
        % read in single experiment
        ind4Exp             = find( strcmp(cribSheet.animal,animal) & cribSheet.experiment_ID == expID );
        expInfo.animal      = cribSheet.animal{ind4Exp(1)};
        expInfo.anNum       = anNum(ind4Exp(1));
        expInfo.Tnames      = cribSheet.trialName(ind4Exp);
        expInfo.fullPath    = fullfile(cribSheet.filePath(ind4Exp),cribSheet.trialName(ind4Exp));
        % additional data (meta data)
        if ~isempty(addFieldNames)
            for j = 1:length(addFieldNames)
                expInfo.(addFieldNames{j}) = cribSheet.(addFieldNames{j})(ind4Exp);
            end
        end
        
        maxNTrial = length(ind4Exp);

    otherwise
        error(['"' method '" is no valid method to read in data. Come on man.' ]);
end

end

