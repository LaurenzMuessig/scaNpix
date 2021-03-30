function index = filterDataGUI(filtVals,criterion,type, direction)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



switch lower(type)
    case 'bylabel'
        index = strcmp(filtVals,criterion);
    case 'spatialinfo'
    case 'custom'
    otherwise
end


end

