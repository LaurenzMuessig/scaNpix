function trialData = batchGetTrialSequence(trialTypes,dirParent)
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
if nargin == 1; dirParent = 'S:\1postDoc\Neuropixels\rawData\'; end

%% 
FolderStruct = dir(dirParent);
trialData    = cell(0,5);

% loop over animal folders
for i = 1:length(FolderStruct)
   
    % . and .. 
    if  ~isfolder(fullfile(FolderStruct(i).folder,FolderStruct(i).name)) || ~isempty(regexp(FolderStruct(i).name,'[.]*','once')); continue; end

    animalFolderStruct = dir(fullfile(FolderStruct(i).folder,FolderStruct(i).name));

    expCounter = 1; % reset counter for current dataset
    % loop over data folders of individual animals

    for j = 1:length(animalFolderStruct)

        if  ~isfolder(fullfile(animalFolderStruct(j).folder,animalFolderStruct(j).name)) || ~isempty(regexp(animalFolderStruct(j).name,'[.]*','once')); continue; end

        tmpData      = getTrialSequence(trialTypes,fullfile(animalFolderStruct(j).folder,animalFolderStruct(j).name));
        tmpData(:,4) = repmat({FolderStruct(i).name},size(tmpData,1),1);   %
        tmpData(:,5) = num2cell(ones(size(tmpData,1),1) .* expCounter,2);  %
        % output
        trialData    = vertcat(trialData,tmpData);
        expCounter   = expCounter + 1; % bump
        %
    end
    %
end
    



end