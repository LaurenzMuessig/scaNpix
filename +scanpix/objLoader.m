function obj = objLoader(objType,dataPath, varargin )
% obj - Load neuropixel or DACQ data from raw into a class object. 
% If you do not want to use defaults you need to generate your own param
% containers (obj params) or struct (map params).
% We only load pos and spike data by default!
% 
%
% Usage:    obj = objLoader( objType, dataPath )
%           obj = objLoader( objType, dataPath, 'inputName', inputVal, .. etc .. )
%
% Inputs:
%    objType   - 'npix' or 'dacq' - data type
%    dataPath  - path to data can be cell array of path strings
%    objParams - optional; containers.Map (see 'scanpix.helpers.defaultParamsContainer' for details on format) or name of file w/o extension (needs to be located in 'Path/To/+scaNpix/files/YourFile.mat')
%    mapParams - optional; mapParamsStruct (see 'scanpix.maps.defaultParamsRateMaps' for details on format) or name of file w/o extension (needs to be located in 'Path/To/+scaNpix/files/YourFile.mat') 
%
% Outputs: obj - class object with data loaded
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TO-DO


%% PARAMS
defaultObjParams = [];
defaultMapParams = [];
loadPos          = true;
loadSpikes       = true;
loadLFP          = false;

p = inputParser;
addParameter(p,'objParams',defaultObjParams, @(x) isa(x,'containers.Map') || ischar(x) || isempty(x));
addParameter(p,'mapParams',defaultMapParams, @(x) isstruct(x) || ischar(x) || isempty(x));
addParameter(p,'pos',loadPos, @islogical);
addParameter(p,'spikes',loadSpikes, @islogical);
addParameter(p,'lfp',loadLFP, @islogical);
parse(p,varargin{:});

loadStr = {'pos','spikes','lfp'};

classFolder = what('+scanpix');

%% Load Data (create class object)
obj = scanpix.ephys(objType,'default',false); % assume default params

% change params if desired
if ~isempty(p.Results.objParams)
    if ischar(p.Results.objParams) 
        if exist(fullfile(classFolder.path,'files',[p.Results.objParams '.mat']),'file') == 2
            obj.changeParams('file', fullfile(classFolder.path,'files',[p.Results.objParams '.mat']) );
        else
            error('scaNpix::objLoader: Can''t find %s .',fullfile(classFolder.path,'files',[p.Results.objParams '.mat']));
        end
    else
        obj.params = p.Results.objParams;
    end
end
% change map params if desired
if ~isempty(p.Results.mapParams)
    if ischar(p.Results.mapParams)
        if exist(fullfile(classFolder.path,'files',[p.Results.mapParams '.mat']),'file') == 2
            tmp = load(fullfile(classFolder.path,'files',[p.Results.mapParams '.mat']));
            f = fieldnames(tmp);
            obj.mapParams = tmp.(f{1});
        else
            error('scaNpix::objLoader: Can''t find %s .',fullfile(classFolder.path,'files',[p.Results.mapParams '.mat']));
        end
    else
        obj.mapParams = p.Results.mapParams;
    end
end

% 
if ischar(dataPath)
    dataPath = {dataPath};
end
% parse directories and trialnames
[filepath,name,~] = cellfun(@fileparts,dataPath,'uni',0);
obj.dataPath = [filepath filesep];

if strcmpi(objType,'npix')
    trialNames = cellfun(@(x) x{1}, cellfun(@(x) regexp(x,'[.]ap','split'),name,'uni',0),'uni',0);
else
    trialNames = name;
end
obj.trialNames = trialNames;

% load
loadStr = loadStr( [p.Results.pos p.Results.spikes p.Results.lfp] );
obj.load(loadStr);



end






