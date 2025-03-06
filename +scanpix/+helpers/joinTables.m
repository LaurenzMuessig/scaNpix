function ResOut = joinTables(ResT1,ResT2,options)
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

%% arguments
arguments
  ResT1 {mustBeA(ResT1,'table')}
  ResT2 {mustBeA(ResT2,'table')}
  options.keys (1,:) {mustBeA(options.keys,'cell')}    = {'dataset','cellID'};
  options.rvars (1,:) {mustBeA(options.rvars,'cell')}  = {};
end
%
if isempty(options.rvars)
    rightVars = ResT2.Properties.VariableNames;
else
    rightVars = options.rvars;
end

%%
[ResOut,iLeft] = innerjoin(ResT1,ResT2,'keys',options.keys,'RightVariables',rightVars);
% keep output in order of left table
[~, sortinds]  = sort(iLeft);
ResOut         = ResOut(sortinds,:);

end