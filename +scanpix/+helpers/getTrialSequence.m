function trialData = getTrialSequence(method,dataDir,trialTypes,options)
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
arguments
  method (1,:) {mustBeMember(method,{'all','any','seq'})} = 'any';
  dataDir (1,:) {mustBeFolder(dataDir)} = 'S:\1postDoc\Neuropixels\rawData\';
  trialTypes (1,:) {mustBeA(trialTypes,'cell')} = {};
  options.exactflag (1,1) {mustBeNumericOrLogical} = true;
  options.excludeflag (1,1) {mustBeNumericOrLogical} = false;
  options.mode (1,:)  {mustBeMember(options.mode,{'pattern','exact'})} = 'pattern';
  options.bslKey (1,:) {mustBeA(options.bslKey,'cell')} = {'fam'};
  options.ignKey (1,:) {mustBeA(options.ignKey,'cell')} = {'sleep'};
  options.getFlankBSL (1,1) {mustBeNumericOrLogical} = true;
end

%%
pathBinFiles = dir(fullfile(dataDir,'**','*.ap.bin'));
pathXMLFiles = dir(fullfile(dataDir,'**','metaData.xml'));
pathXMLFiles = pathXMLFiles(ismember({pathXMLFiles(:).folder},{pathBinFiles(:).folder})); % ignore xml files in session kilosort result folders

if isempty(pathBinFiles) || isempty(pathXMLFiles)
    warning('No raw data found in %s', dataDir);
    trialData = cell(0,4);
    return
end
[~,sortInd]  = sort(datetime({pathBinFiles.date}));
pathXMLFiles = pathXMLFiles(sortInd);
pathBinFiles = pathBinFiles(sortInd);
%%
% need to loop over xml files to get all trial types
currTrialTypes = cell(length(pathXMLFiles),1);
for i = 1:length(pathXMLFiles)
    tmpXML            = scanpix.fxchange.xml_read(fullfile(pathXMLFiles(i).folder,pathXMLFiles(i).name));
    currTrialTypes{i} = lower(tmpXML.trialType);
end

%%
% check if current dataset contains trial sequence
if strcmp(method,'all')
    trialTypes = currTrialTypes;
end
%    
switch method
    case {'any','all'}
        trialInd       = ismember(currTrialTypes, trialTypes);
        if options.excludeflag
            trialInd   = ~trialInd;
            if sum(trialInd) == 0
                trialData = cell(0,4);
                return;
            end
        end
    case 'seq'
        trialInd = scanpix.helpers.matchTrialSeq2Pattern(currTrialTypes,trialTypes,'exactflag',options.exactflag,'bslkey',options.bslKey,'ignKey',options.ignKey,'mode',options.mode,'getFlankBSL',options.getFlankBSL);
        if isempty(trialInd); trialData = cell(0,4); return; end
end

% select relevant trials and sort by date
pathBinFiles   = pathBinFiles(trialInd(1,:));
currTrialTypes = currTrialTypes(trialInd(1,:));
% [~,sortInd]    = sort(datetime({pathBinFiles.date}));
% generate outut data
trialData      = strcat({pathBinFiles(:).folder}, filesep)'; %
trialData(:,2) = {pathBinFiles(:).name};                     %
trialData(:,3) = currTrialTypes;
trialData(:,4) = num2cell(logical(trialInd(2,:)),1);


end

