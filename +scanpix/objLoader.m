function obj = objLoader(objType,dataPath, objParams, mapParams )
% obj - Load neuropixel or DACQ data from raw into a class object. 
% If you do not want to use defaults you need to generate your own param
% containers (obj params) or struct (map params).
% We only load pos and spike data. If you want the LFP data you need to
% load it yourself at later stage with e.g. obj.load('lfp')
% 
%
% Usage:    obj = objLoader( objType, dataPath )
%           obj = objLoader( objType, dataPath, 'inputName', inputVal, .. etc .. )
%
% Inputs:
%    objType   - 'npix' or 'dacq' - class object type
%    dataPath  - path to data can be cell array of path strings
%    objParams - optional; containers.Map (see 'scanpix.helpers.defaultParamsContainer' for details on format) or path to file on disk
%    mapParams - optional; mapParamsStruct (see 'scanpix.maps.defaultParamsRateMaps' for details on format)
%
% Outputs: obj - class object with data loaded
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TO-DO


%% PARAMS

classFolder = what('+scanpix');


%% Load Data (create DACQ object)
if strcmpi(objType,'dacq')
    obj = scanpix.dacq('default',false); % assume default params
elseif strcmpi(objType,'npix')
    obj = scanpix.npix('default',false); % assume default params
else
end
% change params if desired
if nargin > 2 && ~isempty(objParams)
    if ischar(objParams)
        scanpix.helpers.changeParams(obj,'file', fullfile(classFolder.path,'files',[objParams '.mat']) );
    elseif strcmp(classs(objParams),'containers.Map')
        obj.params = objParams;
    end
end
% change map params if desired
if nargin == 4 && ~isempty(mapParams)
    obj.mapParams = mapParams;
end

% 
if ischar(dataPath)
    dataPath = {dataPath};
end
% parse directories and trialnames
[filepath,name,~] = cellfun(@fileparts,dataPath,'uni',0);
obj.dataPath = filepath;

if strcmpi(objType,'npix')
    trialNames = cellfun(@(x) x{1},cellfun(@(x) regexp(x,'[.]','split'),name,'uni',0),'uni',0);
else
    trialNames = name;
end
obj.trialNames = trialNames;

% load
obj.load({'pos','spikes'});



end






