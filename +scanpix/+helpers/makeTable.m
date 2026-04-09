function T = makeTable(varList,varSz,varType,nRows)
% makeTable - preallocate a custom table
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.makeRateMaps(varList, varSz, varType)
%       scanpix.helpers.makeRateMaps(varList, varSz, varType, nRows)
%
% Inputs:
%    varList    - cell array of variable names for table
%    varSz      - column size for each variable name (numeric)
%    varType    - cell array of data type (string) for each variable name (must be
%                 'cell' or 'num')
%    nRows      - n of rows to preallocate (set to 0 if making empty dummy table)
%
% Outputs:
%   T           - preallocated table
%
%
% LM 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
arguments
    varList (:,1) {mustBeA(varList,'cell')}
    varSz (:,1) {mustBeNumeric} 
    varType (:,1) {mustBeA(varType,'cell')}
    nRows (1,1) {mustBeNumeric} = 1;
end


%%
tableData = cell(length(varList),2);
for i = 1:length(varList)
    
    switch varType{i}
        case 'cell'
            scoreAllocated = cell(1,varSz(i));
        case 'num'
            scoreAllocated = nan(1,varSz(i));
        otherwise
    end
    %
    tableData(i,:) = {varList{i},scoreAllocated};  
end

%%
tableData = tableData';
%
T                          = repmat(cell2table(tableData(2,:)),nRows,1); % repmat seems to be the only way to pre-allocate table with n rows and variable column format
T.Properties.VariableNames = tableData(1,:);

end