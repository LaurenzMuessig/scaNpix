function changeParams(obj, changeMode, optionalPath2file)
% changeParams - change params container in obj
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.changeParams
%       scanpix.helpers.changeParams(changeMode)
%       scanpix.helpers.changeParams(changeMode, optionalPath2file)
%
% Inputs:
%    changeMode - 'default' - uses default parameters;
%                 'ui' - opens UI dialogue (default); 
%                 'file' - load params from file
%                 (Note that this gets passed on as 'mode' to scanpix.helpers.getParams)
%
%
%    optionalPath2file - string; optional name of .mat file with parameters (we assume it's in 'YourPathOnDisk\+scanpix\files\')
%                        if ommited and changeMode='file' will open UI dialogue to fetch file
%
% Outputs:
%
% see also: scanpix.helpers.getParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    changeMode = 'ui';
end

if nargin < 3
    optionalPath2file = [];
end

% change params
prmsMap  = scanpix.helpers.getParams(obj,changeMode, optionalPath2file);
% in case we aborted loading we want to keep things as they are
% and don't change anything
if ~isempty(prmsMap)
    obj.params = prmsMap;
end
end
