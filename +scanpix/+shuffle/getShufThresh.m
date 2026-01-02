function [thresh, distr] = getShufThresh(shufVals,options)
% getShufThresh - get percentile from shuffled ditribution
% package: scanpix.shuffle
%
%
% Syntax:
%       [thresh, distr] = scanpix.shuffle.getShufThresh(shufVals)
%       [thresh, distr] = scanpix.shuffle.getShufThresh(shufVals, Name-Value comma separated list)
%
% Inputs:
%    shufVals    - nCellsxnTrials cell array of shuffled values
%    options     - name-value: comma separated list of name-value pairs
%
% Outputs:
%   thresh       - percentile threshold
%   distr        - ditribution of shuffled values 
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 
arguments
    shufVals {mustBeA(shufVals,'cell')}
    options.pctls (1,:) {mustBeNumeric} = 95;
    options.type (1,:) {mustBeMember(options.type,{'cell','pop'})} = 'cell'; 
end

%%
switch options.type
    case 'cell'
        thresh  = cellfun(@(x) prctile(x,options.pctls),shufVals,'UniformOutput',0);
        distr   = [];
    case 'pop'
        distr   = vertcat(shufVals{:});
        thresh  = prctile(distr(:),options.pctls);
end

end