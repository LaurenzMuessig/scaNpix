function filenameOut = checkSaveFile(filenameIn)
% checkSaveFile - check if a file already exists before saving and optionally 
% change the filename 
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
%
% see also:
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

overwriteFile = true; % Default to overwriting or creating new file.
[path,tmpName,ext] = fileparts(filenameIn);

%% check if file exist
if exist(filenameIn, 'file') == 2
  % Ask user if they want to overwrite the file.
  buttonText = questdlg('This file already exists: Do you want to overwrite it?',['Overwrite ''' tmpName ''' ?'], 'Yes', 'No', 'No');
  if strcmpi(buttonText, 'No')
    % User does not want to overwrite. 
    overwriteFile = false;
  end
end
%% edit file name if necessary
if ~overwriteFile
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

