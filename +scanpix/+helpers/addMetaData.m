function addMetaData(obj,name,varargin)
% addMetaData - add meta data tags to object
% These will be added to obj.metaData.(tag)
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.addMetaData
%       scanpix.helpers.addMetaData(name)
%       scanpix.helpers.addMetaData(name, varargin)
%
% Inputs:
%    name     - string/cell array of strings; name(s) of field(s) to be added to obj.metaData. Note that fields in 'obj.defaultMetaDataFields' are initialised upon
%               loading data into object. If ommited will open UI dialogue to fetch data
%    varargin - value/cell array of values; values for fields in name. If single field is added to obj with just 1 single trial, just suply value directly.
%               For other cases (multiple trials and/or fields) format is that varargin{1}{1} == values for name{1}, etc...
%               Note that you need to supply values for each trial in object/field
%
% Outputs:
%
% see also:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InputLengthCheck = true;
if nargin < 2
    % add meta data per UI input
    uiInput = inputdlg(['Metadata field name'; strcat(obj.trialNames,' - value')'],'Please indicate meta data field you want to add',1);
    if isempty(uiInput)
        warning('scaNpix: Cancelled adding meta data to object. It would have been so nice to get some more of it...');
        return;
    end
    name      = uiInput{1};
    values{1} = uiInput(2:end);
    % convert numeric input to numbers - slightly involved
    if all(ismember(cell2mat(values{1}'), '0123456789+-.eEdD'))
        values{1} = num2cell(str2double(string([values{1}])));
    end
end

if ischar(name)
    name = {name};
end

if nargin == 2
    % make field but don't add any values to it (e.g. defaults)
    values = {};
    InputLengthCheck = false;
    
elseif nargin > 2
    if length(name) == 1
        values{:} = varargin{:};
    else
        values = varargin{:};
    end
end

% sanity check input format
if InputLengthCheck && length(name) ~= length(values)
    ME = MException('scaNpix:addMetaData:N_values','Number of values doesn''t match the number of fields to be added to meta data struct');
    throw(ME);
end

% add metadata - first parse input, then assign
for i = 1:length(name)
    if isempty(values)
        tmpVals = {[]};
    else
        tmpVals = values(i);
    end
    
    if iscell(tmpVals{1})
        tmpVals = tmpVals{1};
    elseif length(tmpVals{1}) < length(obj.trialNames)
        tmpVals = repmat(tmpVals,1,length(obj.trialNames));
    else
        tmpVals = num2cell(tmpVals{1});
    end
    % assign to fields
    for j = 1:length(obj.trialNames)
        obj.metaData(1).(name{i}){j} = tmpVals{j};  % assign data here
    end
end
end

