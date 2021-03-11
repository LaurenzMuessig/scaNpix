function [trialNameStrOut, select] = selectTrials(trialNameStrIn,obj)
% selectTrials - select from a list which trials to choose for e.g. loading
% data in to dacq or npix class object
% package: scanpix.helpers
%
%
% Syntax:
%       scanpix.helpers.selectTrials(selectTrials)
%       scanpix.helpers.selectTrials(selectTrials,obj)
%
% Inputs:
%    trialNameStrIn  - char/cell array of filename strings that we want to
%                      choose from
%    obj             - dacq or npix class object (optional) - if supplied
%                      we will remove trial names that are not part of 
%                      selection from trial name list in object  
%
% Outputs:
%    trialNameStrOut - user choice of 'trialNameStrIn'
%
% see also: 
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% skip if only one trial in object
if length(trialNameStrIn) == 1
    return;
end

% UI selection
[select, loadCheck] = listdlg('PromptString','Select which Trial(s) to Include:','ListString',trialNameStrIn,'ListSize',[180 100],'CancelString','Keep All');
if ~loadCheck
    return;
else
    trialNameStrOut = trialNameStrIn(select);
    
    if nargin == 2
        ind            = ismember(obj.trialNames,trialNameStrOut);
        obj.trialNames = obj.trialNames(ind);
        obj.dataPath   = obj.dataPath(ind);
    end
end

end

