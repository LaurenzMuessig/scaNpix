function lfpSpeedFilter = speedFilter2lfp(posSpeedFilter, posFs, lfpFs, nSamplesLFP, options)
% speedFilter2lfp - convert a speed filter (or any logical index) in position samples
% into one in LFP samples. By default the filter gets cleaned up first (fill short gaps
% between running bouts and remove short bouts), as speed often flickers around the
% threshold.
% package: scanpix.lfp
%
% Position sample k covers [(k-1)/posFs, k/posFs) and LFP sample j is at time (j-1)/lfpFs,
% i.e. both start at t = 0 (trial start). Each LFP sample takes the value of the position
% sample it falls into.
%
% Usage:
%       lfpSpeedFilter = scanpix.lfp.speedFilter2lfp(posSpeedFilter, posFs, lfpFs)
%       lfpSpeedFilter = scanpix.lfp.speedFilter2lfp(posSpeedFilter, posFs, lfpFs, nSamplesLFP)
%       lfpSpeedFilter = scanpix.lfp.speedFilter2lfp(posSpeedFilter, posFs, lfpFs, nSamplesLFP, 'maxGapDur', 0.5, 'minBoutDur', 1)
%
% Inputs:   posSpeedFilter - logical vector in position samples (true = valid)
%           posFs          - position sample rate (Hz)
%           lfpFs          - LFP sample rate (Hz)
%           nSamplesLFP    - optional; number of LFP samples. LFP samples beyond the end of
%                            the position data are set to false. Default: length of position
%                            data converted to LFP samples
%
% Options (applied in this order, in position samples, before conversion):
%
%   'maxGapDur'    0.25, - fill gaps (false) between two bouts (true) shorter than this (s); 0 = off
%   'minBoutDur'   0.5,  - remove bouts (true) shorter than this (s); 0 = off
%
%   Use 'maxGapDur', 0, 'minBoutDur', 0 for a plain conversion without clean up
%
% Outputs:  lfpSpeedFilter - logical vector in LFP samples (same orientation as input)
%
% Example:  speedFilt = obj.posData.speed{trialInd} >= 2.5 & obj.posData.speed{trialInd} <= 400;
%           lfpFilt   = scanpix.lfp.speedFilter2lfp(speedFilt, obj.trialMetaData(trialInd).posFs, obj.trialMetaData(trialInd).lfpFs, size(obj.lfpData.lfp{trialInd},2));
%
% See also: scanpix.lfp.getThetaPhase
%
% LM 2026
%%
arguments
    posSpeedFilter {mustBeVector, mustBeNumericOrLogical}
    posFs (1,1) {mustBePositive}
    lfpFs (1,1) {mustBePositive}
    nSamplesLFP (1,1) {mustBeInteger, mustBePositive} = floor(numel(posSpeedFilter) * lfpFs / posFs)
    options.maxGapDur (1,1) {mustBeNonnegative} = 0.25
    options.minBoutDur (1,1) {mustBeNonnegative} = 0.5
end

isRowIn        = isrow(posSpeedFilter);
posSpeedFilter = logical(posSpeedFilter(:));

%% clean up filter
% fill short gaps (only gaps between two bouts, not at start/end of data)
if options.maxGapDur > 0
    [gapStart, gapEnd] = findRuns(~posSpeedFilter);
    fillInd            = gapStart > 1 & gapEnd < numel(posSpeedFilter) & (gapEnd - gapStart + 1) < options.maxGapDur * posFs;
    for i = find(fillInd)'
        posSpeedFilter(gapStart(i):gapEnd(i)) = true;
    end
end
% remove short bouts
if options.minBoutDur > 0
    [boutStart, boutEnd] = findRuns(posSpeedFilter);
    for i = find((boutEnd - boutStart + 1) < options.minBoutDur * posFs)'
        posSpeedFilter(boutStart(i):boutEnd(i)) = false;
    end
end

%% convert to LFP samples
posInd         = floor( (0:nSamplesLFP-1)' .* posFs ./ lfpFs + 1e-9 ) + 1; % pos sample each LFP sample falls into (small tolerance for floating point error at bin edges)
lfpSpeedFilter = false(nSamplesLFP,1);
inPos          = posInd <= numel(posSpeedFilter);
lfpSpeedFilter(inPos) = posSpeedFilter(posInd(inPos));

if isRowIn
    lfpSpeedFilter = lfpSpeedFilter';
end

end

%%
function [runStart, runEnd] = findRuns(x)
% start and end indices of consecutive runs of true in logical column vector x
d        = diff([false; x; false]);
runStart = find(d == 1);
runEnd   = find(d == -1) - 1;
end
