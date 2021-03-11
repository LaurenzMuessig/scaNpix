function [ fNames, validTrodes ] = findDACQFiles( dataDir, identifier, type, removeEXT )
% findDACQFiles - find all DACQ files of certain file type that match some identifier
% See: scanpix.dacqUtils
%
% Syntax:  
%      [ fNames, validTrodes ] = scanpix.dacqUtils.findDACQFiles( dataDir, identifier, type)
%      [ fNames, validTrodes ] = scanpix.dacqUtils.findDACQFiles( dataDir, identifier, type, removeEXT )
%
% Inputs:
%      dataDir     - directory to search
%      identifier  - some identifier string that uniquely identifies dataset (i.e. some part of or whole file name, e.g. part of set file name that will be unique to dataset)
%      type        - type of file (extention) - 'set'; 'pos'; 'tet'; 'eeg'; 'egf'
%      removeEXT   - flag to remove extension from 'fNames' (default=false)
%
% Outputs: 
%      fNames      - (list of) filename(s) that match 'identifier' and 'type'
%      validTrodes - valid tetrodes that were recorded in dataset (only applicable when type='tet')
%
%
% LM 2020


% default - keep extension
if nargin == 3 || strcmp(type,'tet')
    removeEXT = 0;
end

validTrodes = [];

% find specific type
switch type
    case 'set'
        temp = dir(fullfile(dataDir,['*' identifier '*.set']));
        fNames = {temp(:).name};   
    case 'pos'
        temp = dir(fullfile(dataDir,['*' identifier '*.pos']));
        fNames = {temp(:).name};       
    case 'eeg'
        temp = dir(fullfile(dataDir,['*' identifier '*.eeg*']));
        sz   = [temp.bytes];    % Exclude 'dummy' eeg files, which are just a partial text header.
        temp = temp( sz>1000 );  %
        fNames = {temp(:).name};   
    case 'egf'
        temp = dir(fullfile(dataDir,['*' identifier '*.egf*']));
        sz   = [temp.bytes];
        temp = temp( sz>1000 );
        fNames = {temp(:).name};           
    case 'tet'
        % a bit more involved to find all tetrode files
        temp = dir(fullfile(dataDir,['*' identifier '*']));
        allFiles = {temp(:).name};
        tetFileInd = ~cellfun('isempty',regexp(allFiles,'(?<=[.])\d')); % all files ending with an integer
        % output files from KK also have integer after a '.', so we need to remove them from list
        notKKFilesInd = cellfun('isempty',regexp(allFiles,'\w*(clu|fet|klg|fmask)\w*'));
        fNames = allFiles(tetFileInd & notKKFilesInd); % all tetrode file 
        % can also get numeric list of tet numbers 
        if nargout == 2
            splitNames  = regexp(fNames, '[.]', 'split'); % split into name + extension
            splitNames  = vertcat(splitNames{:});
            validTrodes = sort(str2double(splitNames(:,2)))'; % terode numbers that were recorded;
        end  
    otherwise        
        error(['"' type '"' 'is not a valid file type! That''s embarrassing mate'])
end

if length(fNames) == 1
    fNames = fNames{:};
end

% remove file extension, if desired
if removeEXT
    splitNames = regexp(fNames, '[.]', 'split'); % split into name + extension
    splitNames = vertcat(splitNames{:});
    fNames = splitNames(:,1)'; 
end
end
