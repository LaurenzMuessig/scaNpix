function loadMeta(obj, trialIterator, noSyncFlag)
% loadMeta - load meta data for neuropixel data
%
% Syntax:  loadMeta(obj,trialIterator)
%
% Inputs:
%    obj           - ephys class object ('npix')
%    trialIterator - numeric index for trial to be loaded
%
% Outputs:
%
% See also: 
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialIterator (1,1) {mustBeNumeric}
    noSyncFlag (1,1) {mustBeNumericOrLogical} 
end

%%
%
obj.trialMetaData(trialIterator).log.missingSyncsAPStream        = 0;
obj.trialMetaData(trialIterator).log.missingFramesPosStream      = 0;
obj.trialMetaData(trialIterator).log.PosLoadingStats             = NaN(3,2);
obj.trialMetaData(trialIterator).log.frameCountCorruptFromSample = NaN;
obj.trialMetaData(trialIterator).log.nInterpSamplesCorruptFrames = NaN;
obj.trialMetaData(trialIterator).log.SyncMismatchPosAP           = 0;
if isKey(obj.params,'InterpPos2PosFs')
    obj.trialMetaData(trialIterator).log.InterpPos2PosFs         = obj.params('InterpPos2PosFs');
else
    obj.trialMetaData(trialIterator).log.InterpPos2PosFs         = false;
    obj.trialMetaData(trialIterator).log.InterpPosSampleTimes    = {};
end

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
%             obj.trialMetaData(trialIterator).envBorderCoords = [min(obj.trialMetaData(trialIterator).envBorderCoords([1,3])),max(obj.trialMetaData(trialIterator).envBorderCoords([1,3])); min(obj.trialMetaData(trialIterator).envBorderCoords([2,4])),max(obj.trialMetaData(trialIterator).envBorderCoords([2,4]))];
obj.trialMetaData(trialIterator).envBorderCoords = [obj.trialMetaData(trialIterator).envBorderCoords(1:2:end); obj.trialMetaData(trialIterator).envBorderCoords(2:2:end)];

% legacy: older xml files won't contain 'objectPos' field
if ~isfield(metaXMLFile,'objectPos')
    obj.trialMetaData(trialIterator).objectPos = [];
end
%
if ~isfield(metaXMLFile,'posFs')
    obj.trialMetaData(trialIterator).posFs = obj.params('posFs'); % unnescessary really
end

obj.trialMetaData(trialIterator).ppm         = [];
obj.trialMetaData(trialIterator).ppm_org     = [];
% obj.trialMetaData(trialIterator).trackLength = []; % add field to xml?

% load channel map
chanMapInfo = dir(fullfile(obj.dataPath{trialIterator},'*kiloSortChanMap*'));
if isempty(chanMapInfo)
    chanMapName = scanpix.npixUtils.SGLXMetaToCoords_v2(fullfile(obj.dataPath{trialIterator})); % create map on the fly
else
    chanMapName = fullfile(obj.dataPath{trialIterator}, chanMapInfo.name);
end
chanMapStruct = load(chanMapName);
f = fieldnames(chanMapStruct);
for i = 1:length(f)
    obj.chanMap(trialIterator).(f{i})      = chanMapStruct.(f{i});
end
obj.trialMetaData(trialIterator).nChanSort = sum(chanMapStruct.connected);
obj.trialMetaData(trialIterator).nChanTot  = 385;

if isempty(obj.dataSetName)
    if ischar(metaXMLFile.animal)
        ratIDStr    = metaXMLFile.animal;
    else
        ratIDStr    = ['r' num2str(metaXMLFile.animal)];
    end
    %
    anFolderInd     = regexp(obj.dataPath{1},ratIDStr,'once');
    dateStr         = regexp(obj.dataPath{1}(anFolderInd+length(ratIDStr):end),'(?<=(/|\\))\d+(_\d)?(?=(/|\\))','match','once');
    obj.dataSetName = [ratIDStr '_' dateStr];
end

% load sync data
if ~noSyncFlag
    [obj.spikeData(1).sampleT{trialIterator}, obj.trialMetaData(trialIterator).missedSyncPulses] = scanpix.npixUtils.loadSyncData(obj, trialIterator);
    obj.trialMetaData(trialIterator).offSet                                                      = obj.spikeData(1).sampleT{trialIterator}(1); 
    obj.spikeData(1).sampleT{trialIterator}                                                      = obj.spikeData(1).sampleT{trialIterator} - obj.trialMetaData(trialIterator).offSet; 
else
    warning('scaNpix::ephys::loadMeta:Sync data not loaded. I hope you know what you''re doing...');
end
%
end

% spikeGLX meta data %
%             metaDataFile = dir(fullfile(obj.dataPath{trialIterator},'*.ap.meta'));
%             fidMeta  = fopen(fullfile(metaDataFile.folder,metaDataFile.name),'r');
%             C        = textscan(fidMeta, '%[^=] = %[^\r\n]');
%             fclose(fidMeta);
%             obj.trialMetaData(trialIterator).nChanTot = sscanf(C{2}{strcmp(C{1},'nSavedChans')},'%d');
%             obj.trialMetaData(trialIterator).nChanAP  = sscanf(C{2}{strcmp(C{1},'snsApLfSy')},'%d%*%*');
%%%% Do we want to add more info from metafile?? %%%%%%%%%%%%%%%%%

