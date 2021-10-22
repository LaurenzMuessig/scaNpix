function addMetaData(obj,name,values)
% addMetaData - add meta data tags to object. These will be added to 
% obj.trialMetaData.(name)
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.addMetaData
%       scanpix.helpers.addMetaData(name)
%       scanpix.helpers.addMetaData(name, varargin)
%
% Inputs:
%    name     - string, name of field to be added to obj.metaData. If ommited will open UI dialogue to fetch data
%    values   - values for field 'name'. if sing;le value is given for dataset containing multiple trials value will
%               be expanded for all trials.
% Outputs:
%
% see also:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    % add meta data per UI input
    uiInput = inputdlg(['Metadata field name'; strcat(obj.trialNames,' - value')'],'Please indicate meta data field you want to add',1);
    if isempty(uiInput)
        warning('scaNpix: Cancelled adding meta data to object. It would have been so nice to get some more of it...');
        return;
    end
    name   = uiInput{1};
    values = uiInput(2:end);
    % convert numeric input to numbers
    for i = 1:length(values)
        if any(ismember(values{i}, '0123456789+-.eEdD'))
            values{i} = str2num(values{i}); %#okAgrow
        end
    end
end

if nargin == 2
    % initialise empty field 
    values = repmat({[]},1,length(obj.trialNames));
end

if ~iscell(values)
    values = {values};
end

if length(values) < length(obj.trialNames)
    values = repmat(values,1,length(obj.trialNames)); % expand
end

% add metadata
for i = 1:length(obj.trialNames)
   obj.trialMetaData(i).(name) = values{i};
end
end

