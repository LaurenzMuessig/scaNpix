function loadMetaDACQ(obj, trialIterator)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialIterator (1,1) {mustBeNumeric}
end

%%

% first load xml file
metaXMLFileInfo = dir(fullfile(obj.dataPath{trialIterator},['meta*_' obj.trialNames{trialIterator} '.xml']));
%
if isempty(metaXMLFileInfo); return; end
%
metaXMLFile     = scanpix.fxchange.xml_read(fullfile(obj.dataPath{trialIterator},metaXMLFileInfo.name)); %

f = fieldnames(metaXMLFile);
for i = 1:length(f)
    if ~strcmpi(f{i},'comment')
        obj.trialMetaData(trialIterator).(f{i}) = metaXMLFile.(f{i});
    end
end

% obj.trialMetaData(trialIterator).envBorderCoords = [obj.trialMetaData(trialIterator).envBorderCoords(1:2:end); obj.trialMetaData(trialIterator).envBorderCoords(2:2:end)];


end