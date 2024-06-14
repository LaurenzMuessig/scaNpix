function [dataIndex, missingTrials] = matchTrialSeq2Pattern(trialSeqObj,pattern2extract,varargin)  
% matchTrialSeq2Pattern - match a pattern of trial types to the trials run for a specific dataset 
% package: scanpix.helpers
%
% Syntax:  
%    scanpix.helpers.matchTrialSeq2Pattern(trialSeqObj,pattern2extract)
%
% Inputs:
%    trialSeqObj     - cell array of trial type sequence
%    pattern2extract - cell array of trial types to extract 
%    varargin        -  
%
%
% Outputs:
%    dataIndex     - numeric index to extract trials in order of 'pattern2extract'
%    missingTrials - numeric index of missing trial types in 'trialSeqObj'
%
% LM 2024
%
%% Params
bslKey  = 'fam';
%
p = inputParser;
addParameter(p,'bslk',  bslKey, ( @(x) ischar(x) || isstring(x) ));

parse(p,varargin{:});
%%
%
firstNonBslInTrialSeq = find( ~strcmp(trialSeqObj,p.Results.bslk), 1, 'first' );
firstNonBslInPattern  = find( ~strcmp(pattern2extract,p.Results.bslk), 1, 'first' );
%
if isempty(firstNonBslInPattern)
    nBSL = length(pattern2extract);
else
    nBSL = firstNonBslInPattern - 1;
end

if isempty(firstNonBslInTrialSeq)
    dataIndex = find(strcmp(trialSeqObj,p.Results.bslk),nBSL,'last'); % extract the n last baseline trials
    return
end

dataIndex = nan(1,length(pattern2extract));
% Get the n last bsl trials before the first probe.
for i = 1:nBSL
    if firstNonBslInTrialSeq-i < 1
        continue
    else
        dataIndex( nBSL - (i-1) ) = firstNonBslInTrialSeq - i;
    end
end
% now do all trials after the (pre-probe) baseline trials
postBSLTrial2Extr = pattern2extract(firstNonBslInPattern:end);
remainTrials      = trialSeqObj(firstNonBslInTrialSeq:end);
%
c1 = firstNonBslInTrialSeq - 1;
c2 = nBSL + 1;
while ~isempty([postBSLTrial2Extr{:}])
    
    tempInd = find( strcmp( remainTrials,postBSLTrial2Extr{1} ), 1, 'first' );
    if isempty(tempInd)
        postBSLTrial2Extr = postBSLTrial2Extr(2:end);
    else
        dataIndex(c2)     = tempInd + c1;    
        postBSLTrial2Extr = postBSLTrial2Extr(2:end);
        remainTrials      = remainTrials(2:end);
        c1                = c1 + 1;
    end
    c2                    = c2+1;
end
%
missingTrials             = find(isnan(dataIndex));
dataIndex(missingTrials)  = []; 


end