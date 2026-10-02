function [xyOut, lenScale] = rawToFitCoords(obj, trialInd, xy)
% rawToFitCoords - Convert coordinates in raw camera pixels (e.g. object or
% barrier positions from the meta data) into the frame of the processed
% position data (obj.posData.XY)
% package: scanpix.maps
%
% Applies the same steps as were applied to the position data when loading:
%   1) scaling to common ppm (if trialMetaData.PosIsScaled)
%   2) fitting to the environment (if trialMetaData.PosIsFitToEnv{1}), using
%      the transform stored by scanpix.maps.scalePosition in PosIsFitToEnv{3}:
%      XYfit = (XY - origin) .* scale + shift
% For data loaded before the transform was stored, we fall back to only
% subtracting the offset (PosIsFitToEnv{2}), which ignores the rescaling done
% by scalePosition (rect. envs ~1-2%; circular envs can be off substantially).
%
% Syntax:
%       xyOut = scanpix.maps.rawToFitCoords(obj, trialInd, xy)
%       [xyOut, lenScale] = scanpix.maps.rawToFitCoords(obj, trialInd, xy)
%
% Inputs:
%    obj        - ephys class object
%    trialInd   - numeric index of trial
%    xy         - nx2 array of [x y] coordinates in raw camera pixels
%
% Outputs:
%    xyOut      - nx2 array of coordinates in the frame of obj.posData.XY{trialInd}
%    lenScale   - 1x2 [x y] factor to convert lengths (e.g. a radius) from
%                 raw camera pixels to the same frame
%
% see also: scanpix.maps.scalePosition; scanpix.maps.makeOVMap; scanpix.vtc.getBarrCoords
%
% LM 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
    xy (:,2) {mustBeNumeric}
end

%%
metaData = obj.trialMetaData(trialInd);

% scaling to common ppm
if isfield(metaData,'PosIsScaled') && ~isempty(metaData.PosIsScaled) && metaData.PosIsScaled
    scaleFact = metaData.ppm / metaData.ppm_org;
else
    scaleFact = 1;
end
xyOut    = xy .* scaleFact;
lenScale = [scaleFact scaleFact];

% fit to environment
if isfield(metaData,'PosIsFitToEnv') && ~isempty(metaData.PosIsFitToEnv) && metaData.PosIsFitToEnv{1}
    if numel(metaData.PosIsFitToEnv) >= 3 && isstruct(metaData.PosIsFitToEnv{3}) && all(isfield(metaData.PosIsFitToEnv{3},{'origin','scale','shift'}))
        T        = metaData.PosIsFitToEnv{3};
        xyOut    = (xyOut - reshape(T.origin,1,2)) .* reshape(T.scale,1,2) + reshape(T.shift,1,2);
        lenScale = lenScale .* reshape(T.scale,1,2);
    else
        warning('scaNpix::maps::rawToFitCoords:No fit transform stored for trial %i (data loaded with older version?). Only the offset is applied, so coordinates might be slightly off. Reload data to fix this.',trialInd);
        xyOut    = xyOut - reshape(metaData.PosIsFitToEnv{2},1,2);
    end
end

end
