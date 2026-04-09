function overlap = getTraceOverlapScore( BFMask, poPrMask, poPrFRMap )
% Get overlap between barrier field (defined in barrier trial) and post-barrier field
% defined in post-barrier trial, weighted by rate in post-barrier trial.
% Intuitively, we are asking: if there are any new fields in post-barrier trial (compared to
% pre-barrier trial), to what extent to these overlap with the position of new fields created
% by barrier insertion.
%
%      [overlap_score] = bvcTrTraceOverlapScore( BFMask, poPrMask, poPrFRMap )
%
% You can supply the mask inputs as either numerical label arrays, or logical masks, function should
% auto-detect. Rate map should be raw rate, NOT z-scored.

% Only a few lines, but I fiddled with the exact equation a bit, and wanted to maintain consistency
% for different calling funcitons.


% Input check:
if isempty(BFMask) || isempty(poPrMask)
    overlap = nan;
    return
end
if ~islogical(poPrMask)
    poPrMask = poPrMask==1;
end
if ~islogical(BFMask)
    BFMask   = BFMask==1;
end


overlapMask = poPrMask & BFMask;
ovLpPoInBF  = sum( poPrFRMap(overlapMask) ) / sum( poPrFRMap(BFMask) );
ovLPBFInPo  = sum( poPrFRMap(overlapMask) ) / sum( poPrFRMap(poPrMask) );
overlap     = mean( [ovLPBFInPo, ovLpPoInBF] );


end