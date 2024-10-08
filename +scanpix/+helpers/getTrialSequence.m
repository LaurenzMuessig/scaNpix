function trialData = getTrialSequence(trialTypes,dataDir)
% getTrialSequence - extract a trial sequence from a raw data directory
% The idea is that you give a parent directory that contains subfolders
% with the raw data fpr an experimental session. From the xml files we will
% extract what trials in the session match the requested types
%
% Syntax: trialData = getTrialSequence(trialTypes, dataDir)
%
% Inputs:
%    trialTypes - str or cell array of trial identifiers
%    dataDir    - parent folder where raw data is located in
%
% Outputs:
%    trialData  - cell array of datapaths to the raw data that matches the  trial types in the session
%                 {/path/2/raw,bin file name,trial type}
%
% See also: scanpix.helpers.batchGetTrialSequence
%
% LM 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  
pathBinFiles = dir(fullfile(dataDir,'**','*.ap.bin'));
pathXMLFiles = dir(fullfile(dataDir,'**','metaData.xml'));
pathXMLFiles = pathXMLFiles(ismember({pathXMLFiles(:).folder},{pathBinFiles(:).folder})); % ignore xml files in session kilosort result folders

if isempty(pathBinFiles)
    warning('No raw data found in %s', dataDir);
    trialData = cell(0,3);
    return
end
%%
% need to loop over xml files to get all trial types
currTrialTypes = cell(length(pathXMLFiles),1);
for i = 1:length(pathXMLFiles)
    tmpXML            = scanpix.fxchange.xml_read(fullfile(pathXMLFiles(i).folder,pathXMLFiles(i).name));
    currTrialTypes{i} = lower(tmpXML.trialType);
end

%%
% check if current dataset contains trial sequence
if ~all(ismember(trialTypes,currTrialTypes))
    trialData = cell(0,3);
else
    % select relevant trials and sort by date
    trialInd       = ismember(currTrialTypes, trialTypes);
    pathBinFiles   = pathBinFiles(trialInd);
    currTrialTypes = currTrialTypes(trialInd);
    [~,sortInd]    = sort(datetime({pathBinFiles.date}));
    % generate outut data
    trialData      = strcat({pathBinFiles(sortInd).folder}, filesep)'; %
    trialData(:,2) = {pathBinFiles(sortInd).name};                     %
    trialData(:,3) = currTrialTypes(sortInd);
end

end