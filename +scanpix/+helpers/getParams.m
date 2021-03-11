function prmsMap = getParams(obj,mode,path2file)
% getParams - grab params container for dacq or npix class object
% package: scanpix.helpers
%
% Syntax:
%       prmsMap = scanpix.helpers.getParams(obj)
%       prmsMap = scanpix.helpers.getParams(obj, mode)
%       prmsMap = scanpix.helpers.getParams(obj,mode, path2file)
%
% Inputs:
%    obj       - dacq or npix class object
%    mode      - 'default' - uses default parameters (default)
%              - 'ui'      - opens UI dialogue to grab parameters
%              - 'file'    - load params from file
%    path2file - optional name of .mat file with parameters (we assume it's in 'PathOnDisk\+scanpix\files\')
%                if ommited and mode='file' will open UI dialogue to fetch file
%
% Outputs:
%    prmsMap   - params map container
%
% see also: scanpix.helpers.changeParams
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1 || strcmpi(mode,'default')
    
    prmsMap = scanpix.helpers.defaultParamsContainer(class(obj));
    
elseif strcmpi(mode,'ui')
    % gather input from user
    if isempty(obj.params)
        prmsMap = scanpix.helpers.defaultParamsContainer(class(obj)); % we'll use the default params as base and then create a dialogue from these
    else
        prmsMap = obj.params; % we'll use the current ones as base in case user wants to update existing ones
    end
    prompts     = prmsMap.keys;
    defaultVals = prmsMap.values;
    output      = scanpix.helpers.makeCustomUIDialogue(prompts, defaultVals);
    % exit gracefully
    if isempty(output)
        prmsMap = containers.Map;
        return
    end
    % create the container
    prmsMap     = containers.Map(output(:,1),output(:,2));
    
elseif strcmpi(mode,'file')
    % load from file
    if nargin < 3 || isempty(path2file)
        classFolder = what('scanpix');
        [fNames, dataDir] = uigetfile(fullfile(classFolder.path,'files','*.mat'),'Select params containers Map to load.');
        if isnumeric(fNames)
            warning('scaNpix: Loading of params container aborted!');
            prmsMap = containers.Map;
            return;
        end
        path2file         = fullfile(dataDir,fNames);
    end
    temp                  = load(path2file);
    f                     = fieldnames(temp);
    prmsMap               = temp.(f{1});
    
else
    ME = MException('scaNpix:getParams:invalidMode', [mode ' is not an accepted input to construct parameter container. Try ''default'',''ui'' or ''file'' instead.']);
    throw(ME);
end
end

