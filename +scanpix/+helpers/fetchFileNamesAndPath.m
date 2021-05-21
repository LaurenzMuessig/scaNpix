function [trialNames, dataDir] = fetchFileNamesAndPath(fileExt, defaultDir)
% fetchFileNamesAndPath - grab filenames for trials and path on disk for data
% The idea is that we choose a top level directory and then we'll find all
% data files with '.fileExt' in the subdirectories. We will use a UI dialogue
% to fetch that.
%
% Syntax:
%    [trialNames, dataDir] = scanpix.helpers.fetchFileNamesAndPath(fileExt)
%    [trialNames, dataDir] = scanpix.helpers.fetchFileNamesAndPath(fileExt, defaultDir )
%
% Inputs:
%    defaultDir - path to directory where UI dialogue will start
%
% Outputs:
%    trialNames - cell array of data file names (w/o extension)
%    dataDir    - full path(s) to data on disk
%
% see also:
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    classFolder = what('scanpix');
    defaultDir = [classFolder.path filesep];
end

dirList = scanpix.fxchange.uigetfile_n_dir(defaultDir,'Select Top Level Folder(s) With Session Data');
if isempty(dirList)
    warning('scaNpix::fetchFileNamesAndPath:Loading Aborted. Und tschuess!');
    trialNames = [];
    dataDir    = '';
    return
end

% now gather all datasets - search all subfolders for data files
trialStruct = struct(dir('1'));
for i = 1:length(dirList)
    trialStruct = vertcat(trialStruct,dir(fullfile(dirList{i},'**',['*' fileExt])));
end

if isempty(trialStruct)
    warning('scaNpix::fetchFileNamesAndPath:No ''%s'' file(s) found in selcted folder(s). Go and find your data yourself.',fileExt);
    trialNames = [];
    dataDir    = '';
    return
end

[~,ind]   = sort([trialStruct.datenum]);
tempNames = {trialStruct(ind).name}';
dataDir   = strcat({trialStruct(ind).folder}', filesep);

% remove extentions
C = regexp(tempNames,fileExt,'split');
trialNames = cellfun(@(x) x{1}, C,'uni',0);

end


