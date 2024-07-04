function ResOut = joinTables(ResT1,ResT2,varargin)
% joinTables - wrapper for 'innerjoin'. Add data from 'ResT2' to 'ResT1' 
% package: scanpix.helpers
%
%
% Syntax:
%       ResOut = scanpix.helpers.joinTables(ResT1,ResT2)
%       ResOut = scanpix.helpers.joinTables(ResT1,ResT2,key-value pairs)
%
% Inputs:
%       ResT1 - table to add data to
%       ResT2 - table to add data from 
%
% Outputs: 
%
%       ResOut - output table
%
% LM 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Params
keys  = {'dataset','cellID'};
rVars = {};
%
p = inputParser;
addParameter(p,'keys',  keys, ( @(x) iscell(x) || ischar(x) ));
addParameter(p,'rvars', rVars,( @(x) iscell(x) || ischar(x) ));
%
parse(p,varargin{:});

if isempty(p.Results.rvars)
    rightVars = ResT2.Properties.VariableNames;
else
    rightVars = p.Results.rvars;
end

%%
ResOut = innerjoin(ResT1,ResT2,'keys',p.Results.keys,'RightVariables',rightVars);


end