function startGUI(GUIType,classType)
% startGUI - start a GUI to browse various aspects of a data object of
% class dacq or npix
% package: scanpix.GUI
%
% Syntax:
%    scanpix.GUI.startGUI
%
% Inputs:
%    GUIType - string; no argument opens UI dialogue to grab info on what GUI type to start
%            - 'main'       - open main GUI (see scanpix.GUI.mainGUI)
%            - 'lfpBrowser' - open LFP browser GUI (see scanpix.GUI.lfpBrowserGUI)
%            - 'phyHelp'    - open phy helper GUI (see scanpix.GUI.phyHelpGUI)
%
% Outputs:
%
% see also: 
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin == 0
    str = {'main','lfpBrowser'};
    [select, loadCheck] = listdlg('PromptString','Select what GUI to start:','ListString',str,'ListSize',[160 100],'SelectionMode','Single');
    if ~loadCheck
        warning('scaNpix::GUI:startGUI: Seems like you don''t want to use an amazing GUI?');
        return;
    end
    GUIType = str{select};
end

switch lower(GUIType)
    case 'main'
        scanpix.GUI.mainGUI(classType);  
    case 'lfpbrowser'
        % TO DO %
    
    case 'phyhelp'
        [fileN,path] = uigetfile(fullfile(cd ,'*.ap.bin'),'Select a neuropixel data file');
        if isnumeric(path)
            warning('scaNpix::GUI:startGUI: Aborted starting phy helper GUI. That''s lame.?');
            return
        end
        
        npixObj = scanpix.analysis.npix_Loader('dataPath',path,'Tnames',fileN);
        scanpix.maps.addMaps(npixObj,'rate',[],'smooth','boxcar','speedFilterFlagRMaps',1);
        scanpix.maps.addMaps(npixObj,'dir',[]);
        % start UI
        scanpix.GUI.phyHelpGUI(npixObj);
        
    otherwise
        ME = MException('scaNpix:startGUI:InvalidGUIType', ['''' GUIType ''' is not a valid GUI type. You need to write that one yourself I am afraid...' ]);
        throw(ME);
end
end
