function objData = batchLoader(cribSheetPath, method, objType, addMeta, options )
% batchLoader - load a bunch of data based on info from an excel sheet.
% Basically a wrapper for scanpix.objLoader
% See 'scanpix.helper.readExpInfo' for details on sheet format. 
%
% Syntax:
%       objData = scanpix.batchLoader(cribSheetPath, method )
%       objData = scanpix.batchLoader(cribSheetPath, method, Name-Value comma separated list )
%
% Inputs:
%    cribSheetPath - path to cribsheet on disk 
%    method        - 'exp' or 'single'
%    dataType      - 'npix' or 'dacq' - data type
%    varargin      - 'objParams' - containers.Map (see 'scanpix.helpers.defaultParamsContainer' for details on format)
%                  - 'mapParams' - mapParamsStruct (see 'scanpix.maps.defaultParamsRateMaps' for details on format)
%
% Outputs:
%    objData       - cell array with data
%
%
% LM 2021
%
%% Args
arguments
  cribSheetPath (1,:) {mustBeFile}  
  method (1,:) {mustBeMember(method,{'single','singletrial','exp','singleexp'})} = 'exp'
  objType (1,:) {mustBeMember(objType,{'npix','dacq','bhave'})} = 'npix'
  addMeta (1,:) {mustBeA(addMeta,'cell')} = {}
  options.paramsobj (1,:) {mustBeA(options.paramsobj,'containers.Map')} = scanpix.helpers.defaultParamsContainer(objType);
  options.paramsmap (1,:) {mustBeA(options.paramsmap,'struct')} = scanpix.maps.defaultParamsRateMaps;
  options.loadpos (1,1) {mustBeNumericOrLogical} = true;
  options.loadspikes (1,1) {mustBeNumericOrLogical} = true;
  options.loadlfp (1,1) {mustBeNumericOrLogical} = false;
end

% % reformat params
prms      = fieldnames(options)';
prms(2,:) = struct2cell(options);

%%
expInfo = scanpix.helpers.readExpInfo( cribSheetPath, method );

if ~isempty(addMeta) && ~all(ismember(addMeta,fieldnames(expInfo)))
    error('scaNpix::batchLoader: At least one of the fields you want to add as metaData is not included in your experiment info. So this just cannot work...');
end
%
addMetaData = cell(length(addMeta),2);
for i = 1:length(addMeta)
    addMetaData{i,1} = addMeta{i};
    addMetaData{i,2} = expInfo.(addMeta{i});
end
%
objData = cell(length(expInfo), 1);
c = 1;
for i = 1:length(expInfo.animal)

    try
        % load current object
        objData{c} = scanpix.objLoader(objType, expInfo.fullPath{i}, [addMetaData{:,1} cellfun(@(x) x{i},addMetaData(:,2), 'uni',0)], prms{:} );
    catch 
        warning('scaNpix::batchLoader: Couldn''t load dataset from rat %s starting with trial %s',expInfo.animal{i}, expInfo.fullPath{i}{1});
    end
    c = c + 1;
end

end

