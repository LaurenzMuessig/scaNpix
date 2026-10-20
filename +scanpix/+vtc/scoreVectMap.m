function [VMAng,VMDist,VMAngExt,VMDistExt,RP] = scoreVectMap(vMap, options)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here

%%
arguments
    vMap {mustBeNumeric}
    options.binSizeSpat (1,1) {mustBeNumeric} = 2.5;  % cm; spatial bin size of the rate maps the vMap was made from (= vMap distance bin, as vMapBinDist must be 1)
    options.vMapBinAng (1,1) {mustBeNumeric}  = 6;    % deg
end

%%
if all(isnan(vMap))
    [VMAng, VMDist, VMAngExt, VMDistExt] = deal(NaN);
    RP = [];
    return
end

%%
% Shift and wrap the map circularly so that the peak bin is in the centre (deals with issue of field wrapping around map ends).
[~, pkInd]          = max(vMap(:),[],'omitnan');
[~, pkC]            = ind2sub(size(vMap), pkInd);
wrapBin             = pkC+(size(vMap,2)/2);   % wrapBin is the bin that will come to the edge of the wrapped vMap
wrapInd             = mod( (1:size(vMap,2)) + wrapBin, size(vMap,2) );
wrapInd(wrapInd==0) = size(vMap,2);
vMWrap              = vMap(  :,  wrapInd );

% Threshold the map to get main field.
thr = mean(vMap(:),'omitnan') + std(vMap(:),[],'omitnan');

% Get the main field and associated scores.
RP          = regionprops( vMWrap>=thr, vMWrap,  'Convexhull','MaxIntensity', 'WeightedCentroid','BoundingBox' );
[~,pkFdID]  = max(  cat(1,RP.MaxIntensity)  );
RP          = RP( pkFdID );

% 'Unwrap' any scores based on column indices.
RP.ConvexHull(:,1)       = mod(RP.ConvexHull(:,1) + wrapBin, size(vMap,2) );
RP.WeightedCentroid(:,1) = mod(RP.WeightedCentroid(:,1) + wrapBin, size(vMap,2) );

%
pkSubs           = fliplr( RP.WeightedCentroid );
VMDist           = (pkSubs(1)*options.binSizeSpat) - options.binSizeSpat/2;                        % Subtract 1/2 bin to set distance to bin centre not edge.
VMAng            = ( (pkSubs(2)*(options.vMapBinAng*pi/180))-((options.vMapBinAng/2)*pi/180) ) + pi/2;  % radians; subtract 1/2 bin to get bin centre. Orientation of circular axis is E-S-W-N-E
VMAng            = mod(VMAng, 2*pi);  % wrap to [0 2pi)
  %              
VMAngExt         = RP.BoundingBox(3)*(options.vMapBinAng*pi/180);  % field width along angle (radians); widths are unaffected by the wrap
VMDistExt        = RP.BoundingBox(4)*options.binSizeSpat;           % field width along distance


end

% 
% % Get correct distance range
% 
% % Exclude vector bins by dwell.
% if ~isempty( varargin ) && ~isempty( varargin{1} )
%     vMapDwell          = varargin{1};
%     vMap(vMapDwell<0) = nan;
% end