function dataOut = getNExp(ratStr,varargin)
% getNExp - collect the number of exposures to a certain trial type
% 
%
%
% LM 2022
%% TO DO


%% parse varargin
dirParent = 'S:\1postDoc\Neuropixels\rawData\';
trialType = 'fam';
nPre      = 0;


p = inputParser;
addOptional(p,'dirparent',dirParent,@ischar);
addParameter(p,'trialtype',trialType,@ischar);
addParameter(p,'npre',nPre,@isscalar);
parse(p,varargin{:});


%% 
dirIn = fullfile(p.Results.dirparent, ratStr);

FolderStruct = dir(dirIn);
dataOut = cell(0,3);

if isempty(FolderStruct)
    warning('scaNpix::npixUtils::getNExp:: Can''t find a folder called ''%s'' in ''%s''.',ratStr,p.Results.dirparent);
    return
end

c = 1;
for i = 1:length(FolderStruct)
   
    % . and .. 
    if  ~isfolder(fullfile(FolderStruct(i).folder,FolderStruct(i).name)) || ~isempty(regexp(FolderStruct(i).name,'[.]*','once')); continue; end
     
    pathBinFiles = dir(fullfile(dirIn,FolderStruct(i).name,'**','*.ap.bin'));
    pathXMLFiles = dir(fullfile(dirIn,FolderStruct(i).name,'**','metaData.xml'));
    pathXMLFiles = pathXMLFiles(ismember({pathXMLFiles(:).folder},{pathBinFiles(:).folder})); % ignore xml files in session kilosort result folders

    if isempty(pathBinFiles)
        warning('No raw data found in %s%s',FolderStruct(i).folder,FolderStruct(i).name);
        continue
    end
    
    for j = 1:length(pathXMLFiles)
        
        tmpXML = scanpix.fxchange.xml_read(fullfile(pathXMLFiles(j).folder,pathXMLFiles(j).name));
        
        if strcmp(tmpXML.trialType,p.Results.trialtype)
            dataOut{c,1} = pathBinFiles(j).name;
            dataOut{c,2} = pathBinFiles(j).date;
            c = c + 1;
        end
    end
end
% remove empty fields and sort by creation date
indEmpty = cellfun('isempty',dataOut(:,1));
dataOut = dataOut(~indEmpty,:);
[~,sortInd] = sort(datetime(datestr(dataOut(~indEmpty,2))));
dataOut = dataOut(sortInd,:);
dataOut(:,3) = num2cell(1+p.Results.npre:size(dataOut,1)+p.Results.npre,1);
end

