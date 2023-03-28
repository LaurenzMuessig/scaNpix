function [field, value] = grabDefaults_mainGUI(defType)
% grabDefaults_mainGUI - set defaults parameters for plotting etc in main
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
switch lower(defType)
    case 'fig'
        % overview tab
        defStruct.ac_bSz_shortLag  = 0.001;
        defStruct.ac_lag_shortLag  = 0.02;
        defStruct.ac_bSz_longLag   = 0.005;
        defStruct.ac_lag_longLag   = 0.5;
        defStruct.rMapColMap       = 'jet';
        defStruct.rMapColNSteps    = 11;
        defStruct.nWaveForms       = 250;
        defStruct.sACColMap        = 'jet';
        defStruct.OVcolMap         = 'parula';
        defStruct.normColMapByMax  = false;
        defStruct.plotNormSpeedMap = false;
        % compare cells tab
        % linearise tab
    case 'gui'
        defStruct.outputDir        = [cd filesep];
        defStruct.filename_suffix  = datestr(now,'yymmdd');
        defStruct.WF_loadMode      = 'file';
    case 'maps'
        defStruct                     = scanpix.maps.defaultParamsRateMaps;
        defStruct.rate.showWaitBar    = true;
        defStruct.dir.showWaitBar     = true;
        defStruct.objVect.showWaitBar = true;
end
%
field = fieldnames(defStruct);
value = cell(length(field),1);
for i = 1:length(field)
    value{i} = defStruct.(field{i});    
end

end

