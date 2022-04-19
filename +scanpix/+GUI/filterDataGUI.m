function index = filterDataGUI(filtVals, criterion, type, direction)
% filterDataGUI - 
%
% package: scanpix.GUI
%
% Syntax:
%       scanpix.GUI.filterDataGUI(filtVals,criterion,type, direction)
%
% Inputs:
%    filtVals  - 
%    criterion - 
%    type      -
%    direction -           
%
% Outputs:
%    index     - 
%
% see also: 
%
% LM 2020; update 2022



switch lower(type)
    case 'bylabel'
        index = strcmp(filtVals,criterion);
    case 'byvalue'
        fH    = str2func(direction);
        index = fH(filtVals,criterion);
    case 'custom'
    otherwise
end


end

