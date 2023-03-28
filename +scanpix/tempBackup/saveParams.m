function saveParams(obj, type, fileName)
% saveParams - saves current params map container or rate map params to disk as .mat file
% User can save her/his own standard parameters
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.saveParams(type)
%       scanpix.helpers.saveParams(type, fileName)
%
% Inputs:
%    type     - string; 'container' - save general params; 'maps' - save rate map params;
%               Note files should be saved to 'PathOnDisk\+scanpix\files\' and that 'defaultParams' & 'defaultRatemapParams' are
%               not allowed as filenames
%
%    fileName - string; filename for params.mat file
%
% Outputs:
%
% see also:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classFolder = what('scanpix');

switch lower(type)
    case 'container'
        prmsMap = obj.params; %#okagrow
        outStr  = 'prmsMap';
        suggest = 'MyParams';
    case 'maps'
        prms    = obj.mapParams; %#okagrow
        outStr  = 'prms';
        suggest = 'MyRateMapParams';
    otherwise
        ME = MException('scaNpix:saveParams:InvalidType', ['''' type ''' is not a valid type of parameter group you can save. Try ''container'' or ''maps'' instead']);
        throw(ME);
end

if nargin < 3 || isempty(fileName)
    fileName     = {'defaultParams'};
    while any( strcmpi(fileName,{'defaultparams','defaultratemapparams'}) )
        fileName = inputdlg('Please enter filename without extension (''defaultParams'' or ''defaultRateMapParams'' are not allowed as filenames)', 'Save Params', 1, {suggest} );
        if isempty(fileName)
            warning('scaNpix: saving params map container aborted');
            return
        end
    end
else
    if any( strcmpi(fileName,{'defaultparams','defaultratemapparams'}) )
        ME = MException('scaNpix:saveParams:InvalidFilename', '''defaultParams'' or ''defaultRateMapParams'' aren''t allowed as file names mate!');
        throw(ME);
    end
end

save( fullfile(classFolder.path,'files',[ fileName{:} '.mat' ]), outStr);
end

