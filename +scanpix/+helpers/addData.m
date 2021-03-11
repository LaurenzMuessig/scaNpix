function addData(obj,newTrialNames,loadMode,reorderFlag)
% addData - add new data to object, i.e. load additional trials into dataset
% (These will be added to end of sequence in object by default)
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.addData(obj)
%       scanpix.helpers.addData(obj, newTrialNames)
%       scanpix.helpers.addData(obj, newTrialNames, loadMode)
%       scanpix.helpers.addData(obj,newTrialNames, loadMode, reorderFlag )
%       scanpix.helpers.addData([],__)
%
% Inputs:
%    obj           - dacq or npix class object
%    newTrialNames - string(s) of set file to be added to object
%                    if ommited or empty will open UI dialogue to fetch files
%    loadMode      - cell array; {'all'} (default) or any combination of {'pos','spikes','lfp'}
%    reorderFlag   - logical (default=true); if true run also scanpix.helpers.reorderData
%
%
% Outputs:
%
% see also: scanpix.helpers.reorderData
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('scaNpix: All existing rate maps and metadata fields will be removed when adding data to object. You''ll need to re-generate these.'); % could add UI confirmation box?

if ~obj.loadFlag
    warning('scaNpix: You need to load data into the object first before you can use ''scanpix.helpers.addData''.');
    return;
end

if nargin < 2 || isempty(newTrialNames)
    [newTrialNames, newDataDir] = scanpix.helpers.fetchFileNamesAndPath(obj.fileType, obj.params('defaultDir'));
    if isempty(newTrialNames)
        return
    end
    [newTrialNames, ind] = scanpix.helpers.selectTrials(newTrialNames);
end

if nargin < 3
    loadMode = {'all'};
end

if nargin < 4
    reorderFlag = true;
end

if any(strcmp(obj.trialNames,newTrialNames))
    warning('scaNpix: The dataset you want to add is already in object, so it will just be reloaded.');
else
    obj.trialNames = horzcat(obj.trialNames, newTrialNames);
    obj.dataPath   = horzcat(obj.dataPath, newDataDir(ind));
end
obj.load(loadMode, newTrialNames);

% reorder
if reorderFlag
    scanpix.helpers.reorderData(obj);
end

end

