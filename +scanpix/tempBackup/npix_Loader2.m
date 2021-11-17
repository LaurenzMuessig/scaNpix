function NPixObj = npix_Loader( dataPath, objParams, mapParams )
% NPixObj - Load neuropixel data from raw into a npix class object. 
% The parameter space is controlled by 'getParams4NPixLoader' and you can 
% pass any params you want to change as prmsStruct or Name-Value pair list 
% to overwrite the defaults
% Atm we save a params container based on the above to disk and then
% overwrite the defaults after initialising object - this is a bit ugly but
% seems the easiest way
% 
%
% Usage:    NPixObj = npix_Loader
%           NPixObj = npix_Loader( optionalInputStruct )
%           NPixObj = npix_Loader( 'inputName', inputVal, .. etc .. )
%
%
% Outputs:  NPixObj - npix class object with data loaded
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TO-DO


%% PARAMS

classFolder = what('+scanpix');


%% Load Data (create DACQ object)
NPixObj = scanpix.npix('default',false);
if nargin > 1 && ~isempty(objParams)
    if ischar(objParams)
        scanpix.helpers.changeParams(NPixObj,'file', fullfile(classFolder.path,'files',[objParams '.mat']) );
    elseif strcmp(classs(objParams),'containers.Map')
        NPixObj.params = objParams;
    end
end


if nargin == 3
    NPixObj.mapParams = mapParams;
end

% 
if ischar(dataPath)
    dataPath = {dataPath};
end

[filepath,name,~] = cellfun(@fileparts,dataPath,'uni',0);
NPixObj.dataPath = filepath;
trialNames = cellfun(@(x) x{1},cellfun(@(x) regexp(x,'[.]','split'),name,'uni',0),'uni',0);
NPixObj.trialNames = trialNames;

% load
NPixObj.load({'pos','spikes'});



end






