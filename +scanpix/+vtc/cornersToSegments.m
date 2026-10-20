function segs = cornersToSegments(corners)
% cornersToSegments - Convert object corner coordinates (as returned by
% scanpix.vtc.getBarrCoords for 'straight' barriers) into the wall segment
% format needed by scanpix.vtc.makeVectMapBVC
% package: scanpix.vtc
%
% Corners are sorted by angle around their centroid, so the order in which
% they were entered into the meta data doesn't matter (assumes objects are
% convex, e.g. rectangles). With only 2 corners the object is treated as a
% single line (one segment).
%
% Syntax:
%       segs = scanpix.vtc.cornersToSegments(corners)
%
% Inputs:
%    corners - nObj x 2*nPts array, each row [x1..xn y1..yn] (e.g.
%              coords_binned output of scanpix.vtc.getBarrCoords)
%
% Outputs:
%    segs    - nSeg x 5 array, each row [x1 y1 x2 y2 barrID], barrID = row
%              index of object in 'corners' (starts at 1, 0 is reserved for
%              outer walls in makeVectMapBVC)
%
% see also: scanpix.vtc.getBarrCoords; scanpix.vtc.makeVectMapBVC
%
% LM 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
arguments
    corners {mustBeNumeric}
end

%%
if isempty(corners)
    segs = zeros(0,5);
    return
end
if mod(size(corners,2),2) ~= 0
    error('scaNpix::vtc::cornersToSegments:Corner coordinates need to be in format [x1..xn y1..yn] per row.');
end

%%
nPts = size(corners,2)/2;
segs = zeros(0,5);
for i = 1:size(corners,1)
    xy = [corners(i,1:nPts)' corners(i,nPts+1:end)'];
    xy = xy(~any(isnan(xy),2),:);
    if size(xy,1) < 2; continue; end
    % sort corners around centroid so consecutive points are neighbouring corners
    [~, sortInd] = sort( atan2( xy(:,2)-mean(xy(:,2)), xy(:,1)-mean(xy(:,1)) ) );
    xy           = xy(sortInd,:);
    if size(xy,1) > 2
        xyEnd = circshift(xy,-1,1); % close polygon
    else
        xyEnd = xy(2,:);            % line barrier - just one segment
        xy    = xy(1,:);
    end
    segs = [segs; xy xyEnd repmat(i,size(xy,1),1)]; %#ok<AGROW>
end

end
