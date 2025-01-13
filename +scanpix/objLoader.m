function obj = objLoader(objType,dataPath, addMeta, options )
% obj - Load neuropixel or DACQ data from raw into a class object. 
% If you do not want to use defaults you need to generate your own param
% containers (obj params) or struct (map params).
% We only load pos and spike data by default!
% 
%
% Usage:    obj = objLoader( objType, dataPath )
%           obj = objLoader( objType, dataPath, metaData)
%           obj = objLoader( objType, dataPath, metaData, 'inputName', inputVal, .. etc .. )
%           obj = objLoader( objType, dataPath, 'inputName', inputVal, .. etc .. )
%
% Inputs:
%    objType    - 'npix' or 'dacq' - data type
%    dataPath   - path to data can be cell array of path strings
%    metadata   - optional; name by values cell array of metadata to be added to object
%    name-value - comma separated list of name-value pairs   
%                 objParams  - optional; containers.Map (see 'scanpix.helpers.defaultParamsContainer' for details on format) or name of file w/o extension (needs to be located in 'Path/To/+scaNpix/files/YourFile.mat')
%                 mapParams  - optional; mapParamsStruct (see 'scanpix.maps.defaultParamsRateMaps' for details on format) or name of file w/o extension (needs to be located in 'Path/To/+scaNpix/files/YourFile.mat') 
%
% Outputs: obj - class object with data loaded
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TO-DO


%% Args
arguments
  objType (1,:) {mustBeMember(objType,{'npix','dacq','bhave'})}
  dataPath (1,:) {mustBeA(dataPath,'cell')}
  addMeta (:,2) {mustBeA(addMeta,'cell')} = {}
  options.paramsobj (1,:) {mustBeA(options.paramsobj,'containers.Map')} = scanpix.helpers.defaultParamsContainer(objType)
  options.paramsmap (1,:) {mustBeA(options.paramsmap,'struct')} = scanpix.maps.defaultParamsRateMaps
  options.loadpos (1,1) {mustBeNumericOrLogical} = true;
  options.loadspikes (1,1) {mustBeNumericOrLogical} = true;
  options.loadlfp (1,1) {mustBeNumericOrLogical} = false;
end

%%
loadStr = {'pos','spikes','lfp'};
%
% classFolder = what('+scanpix');
%% Load Data (create class object)
obj           = scanpix.ephys(objType,'default',false); % assume default params
obj.params    = options.paramsobj;
obj.mapParams = options.paramsmap;
% 
% parse directories and trialnames
[filepath,name,ext] = cellfun(@fileparts,dataPath,'uni',0);
obj.dataPath        = cellfun(@(x) [x filesep],filepath,'uni',0);
obj.dataPathSort    = cellfun(@(x) [x filesep],filepath,'uni',0);
%
if strcmp(ext,'.dat')
    obj.isConcat = true;
end
%
if strcmpi(objType,'npix')
    trialNames = cellfun(@(x) x{1}, cellfun(@(x) regexp(x,'[.]ap','split'),name,'uni',0),'uni',0);
else
    trialNames = name;
end
obj.trialNames = trialNames;

% load
loadStr = loadStr( [options.loadpos options.loadspikes options.loadlfp] );
obj.load(loadStr);

% add meta data (optional)
if ~isempty(addMeta)
    for i = 1:size(addMeta,1)
        obj.addMetaData(addMeta{i,1},addMeta{i,2});
    end
end



end


% change params if desired
% if isfield(options,'paramsobj')
    % if ischar(p.Results.objParams) 
    %     if exist(fullfile(classFolder.path,'files',[p.Results.objParams '.mat']),'file') == 2
    %         obj.changeParams('file', fullfile(classFolder.path,'files',[p.Results.objParams '.mat']) );
    %     else
    %         error('scaNpix::objLoader: Can''t find %s .',fullfile(classFolder.path,'files',[p.Results.objParams '.mat']));
    %     end
    % else
        % obj.params = options.paramsobj;
    % end
% end
% change map params if desired
% if isfield(options,'paramsmap')
    % if ischar(p.Results.paramsmap)
    %     if exist(fullfile(classFolder.path,'files',[p.Results.mapParams '.mat']),'file') == 2
    %         tmp = load(fullfile(classFolder.path,'files',[p.Results.mapParams '.mat']));
    %         f = fieldnames(tmp);
    %         obj.mapParams = tmp.(f{1});
    %     else
    %         error('scaNpix::objLoader: Can''t find %s .',fullfile(classFolder.path,'files',[p.Results.mapParams '.mat']));
    %     end
    % else
        % obj.mapParams = options.paramsmap;
    % end
% end




