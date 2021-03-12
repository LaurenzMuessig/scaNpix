function startGUI(GUIType,classType)
% startGUI - start a GUI to browse various aspects of the data in obj
% package: scanpix.GUI
%
% Syntax:
%    scanpix.GUI.startGUI
%
% Inputs:
%    GUIType - string; no argument opens UI dialogue to grab info on what GUI type to start
%            - 'main'       - open main GUI (see scanpix.GUI.mainGUI)
%            - 'lfpBrowser' - open LFP browser GUI (see scanpix.GUI.lfpBrowserGUI)
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
        warning('scaNpix: Seems like you don''t want to use an amazing GUI?');
        return;
    end
    GUIType = str{select};
end

switch lower(GUIType)
    case 'main'
        scanpix.GUI.mainGUI(classType);  
    case 'lfpbrowser'
        % TO DO %
        
    %case 'mapgui'
        
    otherwise
        ME = MException('scaNpix:startGUI:InvalidGUIType', ['''' GUIType ''' is not a valid GUI type. You need to write that one yourself I am afraid...' ]);
        throw(ME);
end
end
