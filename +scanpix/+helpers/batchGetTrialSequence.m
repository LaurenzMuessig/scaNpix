function trialData = batchGetTrialSequence(trialTypes,dirParent,options)
% batchGetTrialSequence - extract a trial sequence from a raw data directory
% This is a batch wrapper for scanpix.helpers.getTrialSequence. 'dirParent'
% should be a parent directory that contains folders for individual
% animals. We will then cycle through these and extract all trials that
% the trial types. This is useful to generate a batch loading sheet for a
% specific experimental series.
%
% Syntax: trialData = batchGetTrialSequence(trialTypes, dirParent)
%
% Inputs:
%    trialTypes - str or cell array of trial identifiers
%    dirParent  - parent folder where raw data is located in
%
% Outputs:
%    trialData  - cell array of datapaths to the raw data that matches the  trial types in the session
%                 {/path/2/raw,bin file name,trial type,ratID,experiment group}
%
% See also: scanpix.helpers.getTrialSequence
%
% LM 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO DO 

%%
arguments
  trialTypes (1,:) {mustBeA(trialTypes,'cell')} = {};
  dirParent (1,:) {mustBeFolder} = 'S:\1postDoc\Neuropixels\rawData\';
  options.method (1,:) {mustBeMember(options.method,{'all','any','seq'})} = 'any';
  options.exactflag (1,1) {mustBeNumericOrLogical} = true;
  options.excludeflag (1,1) {mustBeNumericOrLogical} = false;
  options.mode (1,:)  {mustBeMember(options.mode,{'pattern','exact'})} = 'pattern';
  options.bslKey (1,:) {mustBeA(options.bslKey,'cell')} = {'fam'};
  options.ignKey (1,:) {mustBeA(options.ignKey,'cell')} = {'sleep'}; %,'cheeseboard_mem','cheeseboard_of','cheeseboard_task','sqtrack'};
  options.getFlankBSL (1,1) {mustBeNumericOrLogical} = false;
end


%% 
FolderStruct = dir(dirParent);
trialData    = cell(0,6);

% loop over animal folders
for i = 1:length(FolderStruct)
   
    % . and .. 
    if  ~isfolder(fullfile(FolderStruct(i).folder,FolderStruct(i).name)) || ~isempty(regexp(FolderStruct(i).name,'[.]*','once')); continue; end

    animalFolderStruct = dir(fullfile(FolderStruct(i).folder,FolderStruct(i).name));

    expCounter = 1; % reset counter for current dataset
    % loop over data folders of individual animals

    for j = 1:length(animalFolderStruct)

        if  ~isfolder(fullfile(animalFolderStruct(j).folder,animalFolderStruct(j).name)) || ~isempty(regexp(animalFolderStruct(j).name,'[.]*','once')); continue; end

        % tmpData      = scanpix.helpers.getTrialSequence(p.Results.mode,fullfile(animalFolderStruct(j).folder,animalFolderStruct(j).name),trialTypes,'excludeflag',options.excludeflag,'ignoreOrder',options.ignoreOrder,'matchflag',options.matchflag);
        tmpData      = scanpix.helpers.getTrialSequence(options.method,fullfile(animalFolderStruct(j).folder,animalFolderStruct(j).name),trialTypes,'excludeflag',options.excludeflag,'exactflag',options.exactflag,'bslkey',options.bslKey,'ignKey',options.ignKey,'mode',options.mode,'getFlankBSL', options.getFlankBSL);
        tmpData(:,5) = repmat({FolderStruct(i).name},size(tmpData,1),1);   %
        tmpData(:,6) = num2cell(ones(size(tmpData,1),1) .* expCounter,2);  %
        % output
        trialData    = vertcat(trialData,tmpData);
        expCounter   = expCounter + 1; % bump
        %
    end
    %
end
    



end