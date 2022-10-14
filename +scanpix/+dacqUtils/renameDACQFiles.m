function renameDACQFiles
% renameDACQFiles - rename all DACQ files that belong to a specific trial (in case of e.g. typo when saving data)
% We'll grab set file find all other files associated with it and rename them according to edited filename. We'll 
% also make a backup of the original data (which is recommended you do)
%
% LM 2018, 2020 - fixed bug when dealing with output files from running KlustaKwik

%%
makeBackUp = 1;

%%
defaultDir = fileparts(which('renameDACQFiles'));

[fName, dataDir] = uigetfile(fullfile(defaultDir, '*.set'),'Select File From Group That You Want Renamed');
if ~ischar(dataDir)
    warning('Loading was cancelled. Can''t continue. Nooooo.');
    return
end
% get all files from group
fList = dir([dataDir fName(1:end-3) '*']);
% new name
newFName = inputdlg('New file name: ','Please fill in',[1 50],{fName(1:end-4)});
if strcmp(newFName,fName(1:end-4))
     warning('New filename == old filename. This would be a pointless exercise mate, so we better stop here...');
    return
end

% for backup of original data
if makeBackUp
    backupDir = [dataDir 'tempBackup\'];
    if ~isfolder(backupDir)
        mkdir(backupDir);
    end
end
% rename and move
for i = 1:length(fList)
    [~,~,ext] = fileparts(fList(i).name);
    % KlustaKwik files need to be dealt with differently
    if ~isempty( regexp(fList(i).name, '(?<=\.)\w*(?=\.)', 'once') )
        extraInd = strfind(fList(i).name,'.');
        copyfile([dataDir fList(i).name],[dataDir newFName{:} fList(i).name( extraInd(1):extraInd(2)-1 ) ext]);
    else
        copyfile([dataDir fList(i).name],[dataDir newFName{:} ext]);
    end
    movefile([dataDir fList(i).name],[backupDir fList(i).name]);
end

end
