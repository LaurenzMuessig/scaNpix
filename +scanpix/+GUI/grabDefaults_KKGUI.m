function [field, value] = grabDefaults_KKGUI
% grabDefaults_KKGUI - set defaults parameters for plotting etc in main
% GUI. You should edit these if you want to make permanent changes 
% package: scanpix.GUI
%
% Syntax:  
%    scanpix.GUI.grabDefaults_mainGUI(defType)
%
% Inputs:
%    defType - char; indicates which types of default values we want to set
%            - 'fig'  - defaults for plots in GUI window
%            - 'gui'  - general defaults for in GUI
%            - 'maps' - defaults for rate map generation
%
% Outputs:
%    field - cell array with list of field names 
%    value - cell array with corresponding values
%
%
% See also: scanpix.GUI.mainGUI;
%
% LM 2021


% 
defStruct.nWaves      = [];
defStruct.nChans      = 5;
defStruct.nSamples    = 70;
defStruct.precPrePeak = 0.4;

defStruct.outputDir    = '';

%
field = fieldnames(defStruct);
value = cell(length(field),1);

for i = 1:length(field)
    value{i} = defStruct.(field{i});    
end

end
