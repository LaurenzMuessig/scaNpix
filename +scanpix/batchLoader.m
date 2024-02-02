function objData = batchLoader(cribSheetPath, method, dataType, varargin )
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
%%
objParams   = [];
mapParams   = [];
loadPos     = true;
loadSpikes  = true;
loadLFP     = false;
addMetaData = {};     

p = inputParser;
addOptional(p, 'addMetaField', addMetaData, @(x) isempty(x) || iscell(x) || (ischar(x) && ~any(strcmp(x,p.Parameters))));
addParameter(p,'objParams',    objParams,    @(x) isa(x,'containers.Map') || ischar(x) || isempty(x));
addParameter(p,'mapParams',    mapParams,    @(x) isstruct(x) || ischar(x) || isempty(x));
addParameter(p,'pos',          loadPos,      @islogical);
addParameter(p,'spikes',       loadSpikes,   @islogical);
addParameter(p,'lfp',          loadLFP,      @islogical);

parse(p,varargin{:});

% reformat params
prms      = p.Parameters(~strcmp(p.Parameters,'addMetaField'));
tmp       = struct2cell(p.Results);
prms(2,:) = tmp(~strcmp(p.Parameters,'addMetaField'));

%%
expInfo = scanpix.helpers.readExpInfo( cribSheetPath, method );

if ~all(ismember(p.Results.addMetaField,fieldnames(expInfo)))
    error('scaNpix::batchLoader:At least one of the fields you want to add as metaData is not included in your experiment info. So this just cannot work...');
end

if ~iscell(p.Results.addMetaField); addFields = {p.Results.addMetaField}; else; addFields = p.Results.addMetaField; end 
addMetaData = cell(length(addFields),2);
for i = 1:length(addFields)
    addMetaData{i,1} = addFields{i};
    addMetaData{i,2} = expInfo.(addFields{i});
end


objData = cell(length(expInfo), 1);
c = 1;
for i = 1:length(expInfo.animal)

    try
        tmp = cellfun(@(x) x{i}, addMetaData(:,2),'uni',0);
        objData{c} = scanpix.objLoader(dataType, expInfo.fullPath{i}, {addMetaData{:,1}; tmp{:}}', prms{:});
    catch 
        warning('scaNpix::batchLoader:Couldn''t load dataset from rat %s starting with trial %s',expInfo.animal{i}, expInfo.fullPath{i}{1});
    end
    c = c + 1;
end

end

