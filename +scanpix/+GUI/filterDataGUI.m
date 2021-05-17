function index = filterDataGUI(filtVals,criterion,type, direction)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



switch lower(type)
    case 'bylabel'
        index = strcmp(filtVals,criterion);
    case 'byvalue'
        fH = str2func(direction);
        index = fH(filtVals,criterion);
    case 'custom'
    otherwise
end


end

