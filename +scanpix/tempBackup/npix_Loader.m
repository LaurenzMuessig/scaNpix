function NPixObj = npix_Loader( varargin )
% DACQ_dataLoader_v2 - Load DACQ data from raw into a dacq class object. 
% The parameter space is controlled by 'getParams4DACQLoader' and you can 
% pass any params you want to change as prmsStruct or Name-Value pair list 
% to overwrite the defaults
% Atm we save a params container based on the above to disk and then
% overwrite the defaults after initialising object - this is a bit ugly but
% seems the easiest way (could add automatic deletion of file after
% loading?)
% 
%
% Usage:    dacqObj = DACQ_dataLoader
%           dacqObj = DACQ_dataLoader( mode, optionalInputStruct )
%           dacqObj = DACQ_dataLoader( mode, 'inputName', inputVal, .. etc .. )
%
%
% Outputs:  dacqObj - dacq class object with data loaded
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TO-DO

%% PARAMS

prms = scanpix.analysis.getParams4NPixLoader(varargin{:});

% need pos + tet data when making maps
if prms.makeRMaps || prms.makeDirMaps || prms.makeLinMaps 
    prms.loadPos    = true;
    prms.loadSpikes = true;
end

% if prms.loadTet
%     prms.loadPos = true;
% end

classFolder = what('+scanpix');

%% make params container and save to class def folder
keys =  {   'ScalePos2PPM',...
            'posMaxSpeed',...
            'posSmooth',...
            'posHead',...
            'posFs',...
            'loadFromPhy',...
            'APFs',...
            'lfpFs',...
            'defaultDir',...
            'myRateMapParams'};
        
vals =   {  prms.ScalePosPPM,...
            prms.posMaxSpeed,...
            prms.posSmooth,...
            prms.posHead,...
            prms.PosFs,...
            prms.loadFromPhy,...
            3e4,...
            2500,...
            [classFolder.path filesep],...
            prms.path2RateMapParams};
           
prmsMap = containers.Map(keys,vals);

save( fullfile(classFolder.path,'files',[prms.paramsFName '.mat']), 'prmsMap');

%% Load Data (create DACQ object)
switch prms.mode
    case 'hc'
        NPixObj = scanpix.npix('default',false);
        scanpix.helpers.changeParams(NPixObj,'file', fullfile(classFolder.path,'files',[prms.paramsFName '.mat']) ); % ...overwrite with current ones
        NPixObj.dataPath   = prms.dataPath;
        NPixObj.trialNames = prms.Tnames;
    case 'ui'
        NPixObj            = scanpix.npix('ui');    
end
% load
str = {'pos','spikes','lfp'};
str = str([prms.loadPos prms.loadSpikes prms.loadEEG]);
NPixObj.load(str, NPixObj.trialNames);

% % add some useful metadata
% % this is ugly but we want to be flexible so we can add whatever and as
% % many fields as we fancy
% for i = 1:length(prms.metaDataStr)
%     dacqObj.addMetaData(prms.metaDataStr{i},prms.(['val' num2str(i)]));
% end

%% scale path
% WE SHOULD PROBABLY ALLOW INDEXABLE 'minOccForEdge' & 'boxExtent'
if prms.scalePath
    for i = 1:length(prms.triaIndex4Scaling)
        scanpix.maps.scalePosition(NPixObj, prms.triaIndex4Scaling(i), prms.minOccForEdge, prms.boxExtent)
    end
end

%% make maps
% standard rate maps
if prms.makeRMaps
    scanpix.maps.addMaps(NPixObj, 'rate', [], prms );
end

% directional rate maps
if prms.makeDirMaps
    scanpix.maps.addMaps(NPixObj, 'dir', [], prms );
end

% linear rate maps (only when environment is some form of track)
if prms.makeLinMaps
    trackInd = find( ~cellfun('isempty', regexp(lower([NPixObj.trialMetaData.trialType{:}]),'track') ) );
    % NEED TO ADD TRACK PROPS DEFINITION - OR AT LEAST BUILD IN A CHECK
    scanpix.maps.addMaps(NPixObj, 'lin', trackInd, prms ); %#ok                             
end


end






