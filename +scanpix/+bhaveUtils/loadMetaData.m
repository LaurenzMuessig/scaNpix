function loadMetaData(obj, trialIterator)
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
% if ~isfield(metaXMLFile,'objectPos')
%     obj.trialMetaData(trialIterator).objectPos = [];
% end
%
if ~isfield(metaXMLFile,'posFs')
    obj.trialMetaData(trialIterator).posFs = obj.params('posFs'); % HARCODED ATM! Should maybe be added to xml file?
end

obj.trialMetaData(trialIterator).ppm         = [];
obj.trialMetaData(trialIterator).ppm_org     = [];

if isempty(obj.dataSetName)
    if ischar(metaXMLFile.animal)
        obj.dataSetName = [metaXMLFile.animal '_' num2str(metaXMLFile.date)];
    else
        obj.dataSetName = ['r' num2str(metaXMLFile.animal) '_' num2str(metaXMLFile.date)];
    end
end

end

