function filenameOut = checkSaveFile(filenameIn)
% checkSaveFile - check if a file already exists before saving and generate
% an alternative filename in case it does
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.checkSaveFile(filenameIn)
%
% Inputs:
%       filenameIn - 'path/filename.ext'; file to be saved to disk
%
%
% Outputs:
%       filenameOut - alternative output filename to not override original file 
% see also:
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

overwriteFile = true; % Default to overwriting or creating new file.
[path,tmpName,ext] = fileparts(filenameIn);

%% check if file exist and generate unique filename if necessary
if exist(filenameIn, 'file') == 2
  filenameOut = [tmpName '_1'];
    c = 2;
    while exist(fullfile(path,[filenameOut ext]), 'file') == 2
        filenameOut = [tmpName '_' num2str(c)];
        c = c+1;
    end
    filenameOut = fullfile(path,[filenameOut ext]);
else
    filenameOut = filenameIn;
end

end

