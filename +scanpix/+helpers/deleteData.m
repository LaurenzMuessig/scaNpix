function deleteData(obj, trialStr)
% delete - delete data of individual trials from dacq or npix objects.
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.delete
%       scanpix.helpers.delete(trialStr)
%
% Inputs:
%    trialStr - string/cell array of strings; name(s) of trial(s) to be deleted from object
%               if ommited will open UI dialogue to select trial(s)
%
% Outputs:
%
% see also:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse inputs
if nargin < 2
    [select, loadCheck] = listdlg('PromptString','Select what data to delete:','ListString',obj.trialNames,'ListSize',[160 100]);
    if ~loadCheck
        warning('scaNpix: Deleting data aborted. More data is better anyway.');
        return;
    end
    trialStr = obj.trialNames(select);
end

deleteInd = ismember(obj.trialNames,trialStr);
fields = fieldnames(obj);
for i = 1:length(fields)
    if ~any(strcmp(fields{i},obj.fields2spare)) && ~isempty(obj.(fields{i}))
        if isstruct(obj.(fields{i})) && isscalar(obj.(fields{i}))
            obj.(fields{i}) = structfun(@(x) x(~deleteInd),obj.(fields{i}),'uni',0);
        else
            obj.(fields{i}) = obj.(fields{i})(~deleteInd);
        end
    end
end

% remove the last bits if object is now empty
if isempty(obj.trialNames)
    for i = 1:length(obj.fields2spare)
        if ~strcmp(obj.fields2spare{i},'params')
            obj.(obj.fields2spare{i}) = [];
        end
    end
    obj.loadFlag  = false;
    scanpix.helpers.changeParams(obj,'default'); % reset defaults
    obj.mapParams = scanpix.maps.defaultParamsRateMaps(class(obj)); % reset defaults
    warning('scaNpix: Back to square one...');
end
end
