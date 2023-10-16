function startGUI(GUIType,dataType)
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


%
if nargin == 0
    str = {'main','lfpBrowser'};
    [select, loadCheck] = listdlg('PromptString','Select what GUI to start:','ListString',str,'ListSize',[160 100],'SelectionMode','Single');
    if ~loadCheck
        warning('scaNpix::GUI:startGUI: Seems like you don''t want to use an amazing GUI?');
        return;
    end
    GUIType  = str{select};
    dataType = [];
end
%
if nargin == 2
    assert( strcmpi(dataType,'npix') | strcmpi(dataType,'dacq') | strcmpi(dataType,'nexus'), 'scaNpix::GUI:startGUI: ''dataType'' needs to be ''npix'', ''dacq'' or ''nexus''. ' )
end

switch lower(GUIType)
    case 'main'
        scanpix.GUI.mainGUI(dataType);  
    case 'lfpbrowser'
        % TO DO %
    
    case 'phyhelp'
        % start UI
        scanpix.GUI.phyHelpGUI;
        
    otherwise
        ME = MException('scaNpix:startGUI:InvalidGUIType', ['''' GUIType ''' is not a valid GUI type. You need to write that one yourself I am afraid...' ]);
        throw(ME);
end



end
