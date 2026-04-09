function vMap = makeVectMapBVC( posMap, spkMap, barrDataIn, barrType,options)
% Make vector map of firing: firing rate with respect to distance and direction from walls/cues.
%
%           [vMapRate, vMapPos_Spk] = bvcTrVectMap( posMap, spkMap, mapSize, barrStr, barrCirc )
%
% posMap is a position map array, spkMap is a spike map array (or a cell array of spike map, one for each cell recorded in trial).
%
% barrStr format:  [x1 y1, x2 y2,       barrID;   etc   ] defines each straight line of a rectanglar object.
% barrCirc format: [xCen, yCen, radius, barrID;   etc   ] defines each circular object.
%
% 'barrID' is an incrementing integer starting at 1 (1, 2, 3 etc) defining separate physical objects within the arena.
%
%
% vMapRate is a cell array of vector rate maps, one for each recorded cell. 
% vMapPos_Spk  is a cell array of (smoothed) vector position and spike maps, position {1} in array is pos map, {2:nCell+1} are the spike maps.
% 
% Closely based on 'makeBVCMaps_v2' (TW Jan 2018) - uses the same method of defining wall positions with
% respect to binned allocentric position.

%%% IMPORTANT NOTE May 2020: there is a bug in this function, BUT ONLY if prms.vMapBinDist is ~=1.
%%%                          I tried to do a rewrite to allow fine spatial position bins and larger vector
%%%                          bins, but I made a mistake as I forgot that the algorithm relies on a one-to-one correspondence
%%%                          between spatial and vector bins in parts 4a and 4b. BUT works fine when prms.vMapBinDist=1.
%%%                          Need to rethink further binning, could come after 4a-b?
arguments
    posMap {mustBeNumeric}
    spkMap
    barrDataIn {mustBeNumeric}
    barrType (1,:) {mustBeMember(barrType,{'straight','circ'})} = 'straight';
    options.spatBinForVMap (1,1) {mustBeNumeric} = 10;
    options.wallSegOverSamp (1,1) {mustBeNumeric} = 10;   % Relative to spatial map bins.
    options.vMapBinAng (1,1) {mustBeNumeric} = 6;    % Units = deg
    options.vMapBinDist (1,1) {mustBeNumeric} = 1;    % Units = *spatial* map bins. IMPORTANT: see note May 2020 in bvcTrVectMap regarding this parameter (basically must always be 1).
    options.maxEnvSz (1,1) {mustBeNumeric} = 75;  % Largest box dim = 150cm, limit max tuning to half this. ppm=400.
    options.nBinsEdgeThr (1,1) {mustBeNumeric} = 15;
    options.excOuterWall (1,1) {mustBeNumericOrLogical} = 1;
    options.smoothVMap (1,1) {mustBeNumericOrLogical} = 1;
    options.smoothKDist (1,1) {mustBeNumeric} = 5;   % Kernel size for distance, units=vector map bins.
    options.smoothKAng (1,1) {mustBeNumeric} = 30;
end

nBinsDist       = ceil( ( options.maxEnvSz/(0.25*options.spatBinForVMap) )/options.vMapBinDist ) + 1; 
nBinsAng        = 360/options.vMapBinAng;
% ---------------------------------------------------------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0) (New addition). Find the edges of the visited environment, and add these to any barrier co-ords,
%     they should be marked as 'object 0' so following code knows that these are external walls.
if ~iscell(spkMap); spkMap = {spkMap}; end

visMask = posMap>0;
XL      = find( sum(visMask,1)>options.nBinsEdgeThr, 1, 'first' ) - 1;
YL      = find( sum(visMask,2)>options.nBinsEdgeThr, 1, 'first' ) - 1;
if XL == 0 || YL ==0
    posMap                = padarray(posMap,[1 1],NaN,'both');
    spkMap                = cellfun(@(x) padarray(x,[1 1],NaN,'both'), spkMap,'uni',false);
    barrDataIn(:,1:end-1) = barrDataIn(:,1:end-1) + 1;
    visMask               = posMap>0;
    XL                    = find( sum(visMask,1)>options.nBinsEdgeThr, 1, 'first' ) - 1;
    YL                    = find( sum(visMask,2)>options.nBinsEdgeThr, 1, 'first' ) - 1;
end
%
XU      = find( sum(visMask,1)>options.nBinsEdgeThr, 1, 'last'  ) + 1;
YU      = find( sum(visMask,2)>options.nBinsEdgeThr, 1, 'last'  ) + 1;
%
% barrData = [XL,YL, XL,YU, 0; ...
%            XL,YU, XU,YU, 0; ...
%            XU,YU, XU,YL, 0; ...
%            XU,YL, XL,YL, 0; ...
%            barrDataIn];
barrData = [XL,YL, XL,YU, 0; ...
           XL,YU, XU,YU, 0; ...
           XU,YU, XU,YL, 0; ...
           XU,YL, XL,YL, 0];
%
mapSize  = size(posMap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Make the wall and barrier XY coord list. (convert corner coords or centre/radius specs to lines).
%     The main output from this section is 'barrCoords', with dimensions ( nWallPoints, 4 ), where the format of each line is:
%
%           [x, y, barrID, wallID]
%
%     'barrID', is supplied by user, corresponds to individual physical barrier or object. barrID=0 is 
%     reserved for the outer wall.
%     'wallID', is 0 for all outer walls, and then increments 1 per actual wall segment otherwise.
%     Used for removing non-visible points. (Note that this system assumes that outer walls never block view of themselves, though).
barrCoords  = ones( 0, 4 );
wallIDCount = 1;   % Count for unique ID of wall segments within env (i.e. excluding outer wall)
wallIDType  = {};  % Keeps a record of whether internal wall segments of straight or circle.

switch barrType
    case 'straight'
        barrData = [barrData; barrDataIn];
        % 1a. Straight wall segments
        for ii=1:size( barrData,1 )
            barrLength = sqrt( diff( barrData(ii,[1 3]) ).^2 + diff( barrData(ii,[2 4]) ).^2  );  % First work out length then assign number of points of line,
            nPoints    = round( barrLength*options.wallSegOverSamp );                                                % this way, keep even 1/10 grid spacing, even for diagonals.
            barrXTemp  = linspace( barrData(ii,1), barrData(ii,3), nPoints )';
            barrYTemp  = linspace( barrData(ii,2), barrData(ii,4), nPoints )';
            if barrData(ii,5)==0
                wallIDTemp = zeros(size(barrXTemp));
            else
                wallIDTemp              = ones(size(barrXTemp)) .* wallIDCount;
                wallIDType{wallIDCount} = 'S';
                wallIDCount             = wallIDCount + 1;
            end
            barrCoords     = cat(1, barrCoords, [barrXTemp, barrYTemp, ones(size(barrXTemp)).*barrData(ii,5), wallIDTemp]);
        end
    case 'circ'
        % 1b. Circle wall segments.
        %     TODO - need to control that user input is OK, that there is an outer wall (or it fails gracefully if not), that if there are barrStr and barrCirc
        %     inputs, the barrIDs are not conflicting or overlapping.
        for ii=1:size( barrDataIn, 1 )
            nPoints = round( 2*pi*barrDataIn(ii,3)*options.wallSegOverSamp );  % nPoints is circumference (in bins) x10.
            [barrXTemp, barrYTemp] = pol2cart( linspace( 2*pi/nPoints, 2*pi, nPoints )', ones(nPoints,1).*barrDataIn(ii,3) );  % First angular binin circle is not 0, so 0/2pi are not double counted.
            barrXTemp              = barrXTemp + barrDataIn( ii, 1 );
            barrYTemp              = barrYTemp + barrDataIn( ii, 2 );
            if barrDataIn(ii,4)==0
                wallIDTemp = zeros(size(barrXTemp));
            else
                wallIDTemp              = ones(size(barrXTemp)) .* wallIDCount;
                wallIDType{wallIDCount} = 'C';
                wallIDCount             = wallIDCount + 1;
            end
            barrCoords     = cat(1, barrCoords, [barrXTemp, barrYTemp, ones(size(barrXTemp)).*barrData(ii,4), wallIDTemp]);
        end
end

nBarrWalls = wallIDCount - 1;
% 1c. If wall segments are continuous, there will be an overlapping point where they meet: remove these by getting uniqueXY rows.
[barrXYUni, ind] = unique( barrCoords, 'rows', 'stable' );  % doinf unique on barrCoords without a row index only gets the corners withon a wall ID, i.e. for the outer wall.
barrCoords       = [barrXYUni, barrCoords(ind)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Some more prep work on the map: (a) work out which map bins the barriers/walls fall in, 
%     so as we can exclude these bins from the analysis, (b) find the convex hull of the outer 
%     wall (that with barrID assigned to 0), and exclude bins outside of this from the analysis.
visMap = true( mapSize );   % 'visMap' will be a logical mask of the bins in the map for which we will create a BVC response.
% 2a. Find barrier bins in map
barrBinsInd          = sub2ind( mapSize, round(barrCoords(:,2)), round(barrCoords(:,1)) );  % Bin<->coord convention is that integers represent bin centres, so 'round' is the correct operation. 
visMap( barrBinsInd) = false;
% 2b. Exclude bins outside outer wall (a procedure which should be avoidable for square, but is essential for circle), 
%     and also those that are inside barriers of a >1bin thickness.
barrIDList = unique( barrCoords(:,3) );
[X,Y]      = meshgrid( 1:mapSize(2), 1:mapSize(1) );
for ii=1:length( barrIDList )
    indBarr    = barrCoords(:,3) == barrIDList(ii);
    coordsTemp = barrCoords( indBarr, : );
    indCH   = convhull( coordsTemp(:,1), coordsTemp(:,2) );
    inMask  = inpolygon( X, Y, coordsTemp(indCH,1), coordsTemp(indCH,2) );
    if barrIDList(ii)==0
        visMap(~inMask) = false;  % If the barrier in question is the outer wall, then set the visMap *outside* to false.
    else
        visMap(inMask) = false;   % If it is an internal barrier, then set the visMap *inside* the barrier to false.
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) For each spatial bin, transform spatial cartesian wall coordinates into bin-centred, 
%     allocentric polar coordinates. The key variable generated is an array showing logical presence/absence
%     of wall, at each dist/ang combination, for each visited cartesian spatial bin ('wallPos_RxThxSB').
%
%     Why is there so much code for such a simple transform? Mainly due to checking for occulsions,
%     so we only get the closest wall when there are multiple walls/barriers (sections 3b and 3c).
%     Should think of a more elegant way of doing this.

[visY,  visX]   = ind2sub( mapSize, find(visMap) );
wallPos_RxThxSB = zeros( nBinsDist, nBinsAng, length(visX) );
for itSB=1:length(visX)  % itSB=iterator over visited Spatial Bins
    
    % 3a. Get a circularly sorted list of barrier position coords (polar), *relative to visited bin itSB* %
    [Th, R]           = cart2pol( barrCoords(:,1)-visX(itSB), barrCoords(:,2)-visY(itSB) );
    [ThSort, sortInd] = sort( Th, 'ascend' );
    RSort             = R( sortInd );
    wallIDSort        = barrCoords( sortInd, 4 );
    
    % 3b. If we are looking stright down a wall from the end (i.e. wall is self occluding), only take the closest segment.
    % Code here assumes that there is only one self-occluding wall at any one time, otherwise won't work.
    diffTh                = diff( ThSort );
    if sum( diffTh==0 ) > 10      % There might always be a few pairs of identical thetas (from different walls), but continuous self-occluding wall must have more than 10. 
        constThInd         = diffTh==0 & diff(wallIDSort)==0;  % The second test here is to make sure it is continuous theta from a continuous wall.
        RSortTemp          = RSort;
        RSortTemp( ~constThInd ) = nan;
        [~, closestSegInd] = nanmin( RSortTemp );
        constThInd(  closestSegInd  ) = false;
        ThSort         = ThSort( ~constThInd );
        RSort          = RSort( ~constThInd );
        wallIDSort     = wallIDSort( ~constThInd );
    end
    
    % 3bii. Deal with circular objects in mid-arena. Need to remove all edge points on the occluded side of the circle.
    %       To acheive this find 'distance to horizon' (where current spatial bin is observer) using secant-tangent theorem
    %       and delete wall points further away than that.
    for itWl=1:nBarrWalls
        if strcmp(wallIDType{itWl},'C')
            O_F            = max( RSort(wallIDSort==itWl) );
            O_C            = min( RSort(wallIDSort==itWl) );
            O_H            = sqrt( O_F.*O_C );
            obscInd        = wallIDSort==itWl & RSort>O_H;
            ThSort         = ThSort( ~obscInd );
            RSort          = RSort( ~obscInd );
            wallIDSort     = wallIDSort( ~obscInd );
        end 
    end
    
    % 3c. When there is a barrier in the middle of the box, need a routine to get the *closest* wall at each angle %
    %     The approch is based on the assumption that walls forwhich the closest point to the rat is the least occulde
    %     others which have the closest point further away. Ths is fine for a simple rectangle barrier, but will not
    %     work for many situations with more complex shapes, or with two barriers.
    % First, get the closest distance to the rat for all of the internal (barrier) walls. %
    closestWallDists = nan(1,nBarrWalls);
    for itWl=1:nBarrWalls
        closestWallDists(itWl) = min( RSort( wallIDSort==itWl ) );
    end
    [closestWallDistsSorted,closestWallOrder]      = sort(closestWallDists,'ascend');
    % Debug stop %
    if 0 && visY(itSB)==3 && visX(itSB)==14
        disp('dum');
    end
    % Then, starting from the closest wall, remove points from all other walls that lie within the extent of that wall %
    for itWl=1:length(closestWallOrder)
        wallNumInd   = find( wallIDSort==closestWallOrder(itWl) );
        if length(wallNumInd) <= 1
            continue  % In the case that wall N has already been completely removed, or having just one segment left (which cant occlude anything).
        end        
        % The only tricky bit here comes from deciding if the wall is lying across the angular wrap point (as the wall segments
        % will never be a completely coherent block, they will always be interspersed with segments of other occluding or occluding wall). 
        wallGaps      = diff( wallNumInd );
        wallExtentInd = false( size(wallIDSort) );
        gapTol        = pi*(3/2);
        if max(wallGaps)<gapTol
            wallExtentInd( wallNumInd(1) : wallNumInd(end) ) = true;    % In the case of a coherent chunk in the middle of the circular FOV ..
        else
            wallExtentInd( 1 : wallNumInd(    find(wallGaps>gapTol,1,'first')    ) )  = true;  % If the wall lies across the wrap point
            wallExtentInd( wallNumInd(   find(wallGaps>gapTol,1,'first')+1   ) : end) = true;
        end
        % Now remove wall segments not from wall N, lying within the extent of wall N %
        segRemoveInd = wallExtentInd & wallIDSort~=closestWallOrder(itWl);
        % Check if the closest point is on a corner of two walls: if it is, don't remove points from the other wall that makes up the corner, 
        % so as to preserve the corner (otherwise, if the walls have sorted the 'wrong' way, we can remove the actual corner without meaning to.
        if itWl==1 && length(itWl)>1 && closestWallDistsSorted(itWl)==closestWallDistsSorted(itWl+1)
            segRemoveInd( wallIDSort==closestWallOrder(itWl+1) ) = false;
        end
        % Now actually remove the segments.
        ThSort       = ThSort( ~segRemoveInd );
        RSort        = RSort( ~segRemoveInd );
        wallIDSort   = wallIDSort( ~segRemoveInd );
    end
    % If requested, remove all points related to the outer wall.
    % This is to allow vector tuning to a discrete onbject ONLY to be measured. 
    % Do this so late in code in case there are other places where existence of an outer wall is assumed.
    if options.excOuterWall
        RSort( wallIDSort==0 ) = nan;   % Set to NaN, so calcualtion of remaining distances can continue unchanged.
    end
    
    % 3d. Now, get the distribution of vectors to wall/barrier allocentric relative to rat, *at this bin XY*.
    %     Bin the wall angles (into final vector map bins), and get mean dist to closest wall/barrier at each angle.
    ThSort(ThSort<=0)    = ThSort(ThSort<=0) + (pi*2);   % Range thus far -pi:pi, shift to 0:2pi
    ThSort               = (ThSort./pi) .* 180;
    ThBin                = ceil( ThSort./options.vMapBinAng );  
    meanRInThBin         = accumarray( ThBin, RSort, [nBinsAng 1], @nanmean, nan );   
    
    % 3e. Turn this into a 2D 'logical histogram' of wall existence (dist x dir, values are 0 or 1 [wall present or not).
    binRInThBin                                               = ceil( meanRInThBin./options.vMapBinDist );   % Bin the distances.
    ThBinExcInd                                               = binRInThBin>nBinsDist | isnan(binRInThBin);
    binRInThBin(ThBinExcInd)                                  = nBinsDist+1;                     % To deal with 'out of range' walls and angle bins where dist is NaN (no walls),
    angDistHistDum                                            = zeros( nBinsDist+1, nBinsAng );  % set these to dist limit +1 .. 
    angDistHistDum( sub2ind( size(angDistHistDum), binRInThBin, (1:nBinsAng)' ) ) = 1;
    wallPos_RxThxSB( :, :, itSB )                             = angDistHistDum( 1:end-1, : );    % .. but remove these 'out of range' counts before final assignment.
end

% % Debug - plot wall pos for each spatial bin.
% if 0
%     figure;
%     if 0
%         for itSB=1:150
%             ax = subplot(10,15,itSB);
%             imagesc(wallPos_RxThxSB(:,:,itSB+500));
%             axis(ax,'off');
%         end
%     else
%         for itDB=1:60
%             ax = subplot(6,10,itDB);
%             imagesc(squeeze(wallPos_RxThxSB(:,itDB,:)));
%             axis(ax,'off');
%         end
%     end
% end
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Generate polar wall dwell and rate maps. Use the list of distances to wall in each allocentric 
%     angular bin, at each cartesian spatial bin as a look up table, to get overall mean based
%     on rat's actual dwell.
% 4a. Dwell map. Linearise unsmoothed pos map, multiply presence/absence wall array by this, so that
%     overall polar map of dist to wall is consistent with rat's actual dwell.
linPosMap    = posMap( visMap );
wallHistWtd  = bsxfun( @times, wallPos_RxThxSB, shiftdim(linPosMap,-2) );
vMap{1}      = sum( wallHistWtd, 3, 'omitnan' );
visMask      = vMap{1}>0;
% 4b. Spk map. Same as above.
%     IMPORTANT: I have played around with this function to get 'dwell maps' from the 'pos map' output, and 'spike maps' that
%                are based solely on rate maps. Upshot is a mismatch between summing 'wallHistWtd' in 4a and meaning 'wallHistWtd'
%                in 4b. This means you *can't* use the function as is for creating rate maps from 'raw' pos and spikes.
for itCl = 1:length(spkMap)
    if sum( spkMap{itCl}(:),'omitnan' ) > 0 % Need to protect against empty inputs: pre-probe only cells in normal datasets.
        linSpkMap    = spkMap{itCl}( visMap );
        wallHistWtd  = bsxfun( @times, wallPos_RxThxSB, shiftdim(linSpkMap,-2) );
        vMap{itCl+1} = mean( wallHistWtd, 3, 'omitnan' );
    else
        vMap{itCl+1} = nan(size(vMap{1}));
    end
end
% 4c. Smooth.
if options.smoothVMap
    % First make filter: guassian,different size+sigma for distance and dir, sigma always 1/3 kernel length.
    fHori         = fspecial( 'gaussian', [1    options.smoothKAng/options.vMapBinAng], (options.smoothKAng/3)/options.vMapBinAng );
    fVert         = fspecial( 'gaussian', [ceil(options.smoothKDist/options.vMapBinDist)  1], (options.smoothKDist/3)/options.vMapBinDist );
    k             = fVert*fHori;
    % Smooth
    visMaskSm = imfilter( repmat(double(visMask),1,3), k );  % Pad for circular smooth of angle (see below).
    for ii=1:length(vMap)
        if ~all(isnan(vMap{ii}))
            padIm    = repmat( vMap{ii}, 1, 3 );
            padIm    = imfilter( padIm, k ) ./ visMaskSm;
            vMap{ii} = padIm(    :,    (1:size(vMap{ii},2)) + size(vMap{ii},2)   );
        else
            vMap{ii} = nan(size(vMap{1}));
        end
    end
end
for ii=1:length(vMap);   vMap{ii}( ~visMask ) = NaN;   end

end
% 4d. Rate maps.
% 27/03/2020: currently not using these.
% vMapR = cell(size(vMap));
% for ii=1:(length(vMap)-1)
%     vMapR{ii} = vMap{ii+1} ./ vMap{1};
% end