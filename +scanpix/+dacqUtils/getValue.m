function rtnVal = getValue(txt,keyStr)
% getValue - Get a value from value-key cell array
% See: scanpix.dacqUtils
%
% Syntax: 
%    rtnVal = scanpix.dacqUtils.getValue(txt,keyStr)
%
% Inputs:
%    txt    - cell array with 'name' (col1) - value (col2) pair 
%    keyStr - 'name' whose value you want to extract 
%
% Outputs:
%    rtnVal - value corresponding to 'name' 
%
% TW (org. Scan function)

ind = strcmp(keyStr,txt(:,1));
if ~any(ind)
    rtnVal = [];
else
    rtnVal = txt{ind,2};
end

end

