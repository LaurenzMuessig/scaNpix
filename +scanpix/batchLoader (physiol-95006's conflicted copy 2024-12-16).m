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
  options.paramsobj (1,:) {mustBeA(options.paramsobj,'containers.Map')} = scanpix.helpers.defaultParamsContainer(objType)
  options.paramsmap (1,:) {mustBeA(options.paramsmap,'struct')} = scanpix.maps.defaultParamsRateMaps
  options.loadpos (1,1) {mustBeNumericOrLogical} = true;
  options.loadspikes (1,1) {mustBeNumericOrLogical} = true;
  options.loadlfp (1,1) {mustBeNumericOrLogical} = false;
end


% objParams   = [];
% mapParams   = [];
% loadPos     = true;
% loadSpikes  = true;
% loadLFP     = false;
% addMetaData = {};     
% 
% p = inputParser;
% addOptional(p, 'addMetaField', addMetaData, @(x) isempty(x) || iscell(x) || (ischar(x) && ~any(strcmp(x,p.Parameters))));
% addParameter(p,'objParams',    objParams,    @(x) isa(x,'containers.Map') || ischar(x) || isempty(x));
% addParameter(p,'mapParams',    mapParams,    @(x) isstruct(x) || ischar(x) || isempty(x));
% addParameter(p,'pos',          loadPos,      @islogical);
% addParameter(p,'spikes',       loadSpikes,   @islogical);
% addParameter(p,'lfp',          loadLFP,      @islogical);
% 
% parse(p,varargin{:});
% 
% % reformat params
prms      = fieldnames(options)';
prms(2,:) = struct2cell(options);

%%
expInfo = scanpix.helpers.readExpInfo( cribSheetPath, method );

if ~isempty(addMeta) && ~all(ismember(addMeta,fieldnames(expInfo)))
    error('scaNpix::batchLoader:At least one of the fields you want to add as metaData is not included in your experiment info. So this just cannot work...');
end

% if ~iscell(p.Results.addMetaField); addFields = {p.Results.addMetaField}; else; addFields = p.Results.addMetaField; end 
addMetaData = cell(length(addMeta),2);
for i = 1:length(addMeta)
    addMetaData{i,1} = addMeta{i}:)
    addMetaData{i,2} = expInfo.(addMeta{i});
end


objData = cell(length(expInfo), 1);
c = 1;
for i = 1:length(expInfo.animal)

    try
        if length(addMetaData{i,2}) == 1 &&  length(expInfo.fullPath{i}) > 1
            addMetaData{i,2} = {repmat(addMetaData{i,2},1,length(expInfo.fullPath{i}))};
        end
        tmp = cellfun(@(x) x{i}, addMetaData(:,2),'uni',0);
        objData{c} = scanpix.objLoader(objType, expInfo.fullPath{i}, {addMetaData{:,1}; tmp{:}}', prms{:} );
    catch 
        warning('scaNpix::batchLoader:Couldn''t load dataset from rat %s starting with trial %s',expInfo.animal{i}, expInfo.fullPath{i}{1});
    end
    c = c + 1;
end

end

