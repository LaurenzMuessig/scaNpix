function [linPos, dirInd, cohRunInd] = linearisePosData(obj, trialInd)
% linearisePosData - linearise postion data from linear track recordings.
% Atm square track or classic linear track is supported.
% package: scanpix.maps
%
%
% Usage:
%       [linPos, dirInd, cohRunInd] = scanpix.maps.linearisePosData(xy,direction,trackType);
%       [linPos, dirInd, cohRunInd] = scanpix.maps.linearisePosData(xy,direction,trackType, prmsStruct);
%       [linPos, dirInd, cohRunInd] = scanpix.maps.linearisePosData(xy,direction,trackType, Name-Value comma separated list);
%
%
% Inputs:   xy          - array of xy positions
%           direction   - array of head directions
%           trackProps  - struct with track propereties. needs to include fields 'type' ('sqtrack' or 'lintrack'), 
%                         'length' (length in cm; for sq track length of 1 arm), 'ppm' (pix/m from tracking)
%                         and 'posFs' (position sample rate)
%           varargin    - prmsStruct: structure with parameter fields to be changed from defaults
%                       - name-value: comma separated list of name-value pairs     
%
% Outputs:  linPos    - lineraised positions
%           dirInd    - index to separate runs in different directions (1=CW and 2=CCW)
%           cohRunInd - logical index to flag coherent runs in data (definded by prms.durThrCohRun)
%
% LM 2020
% 

%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (:,1) {mustBeNumeric}
end

%%
trackLength = obj.trialMetaData(trialInd).ppm * max(obj.trialMetaData(trialInd).envSize)/100; % in pix
% data from object
xy          = obj.posData.XY{trialInd};
HDirections = obj.posData.direction{trialInd};

% positions(addPosFilter,:) = NaN;
%%
if contains(lower(obj.trialMetaData(trialInd).trialType),'sq')
    % 1) Linearise the positions %
    % 1a. First need a good estimate of env edges, for radial-isation routine to work well. Use the same as for box scaling,
    %     i.e. the first camera pixel in each dimension with >=1 sec of summed dwell.
    envEdges = nan(2,2);
    for i = 1:2
        d              = xy(:,i);
        dHist          = histcounts( d, 0:1:ceil(max(d)) );
        goodD          = dHist >= obj.mapParams.linpos.minDwellForEdge * obj.trialMetaData(trialInd).posFs;
        envEdges(:,i)  = [find(goodD,1,'first'); find(goodD,1,'last')];
    end
    % 1b. Radialise square.
    [linPos, armDims]  = radialLineariseSquare_v2(xy,envEdges,trackLength);
    linPos             = linPos + 1;  % Convert from 0-based to 1-based
    %     armDims            = armDims + 1; % for ease of indexing etc later on.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2 Sort out the CW and CCW runs properly ('radialLineariseSquare_v2' only does it on a sample-by-sample basis) %
    % First, rotate the dirs 90deg CCW, so as to switch from phys/eng convention to DACQ convention - makes the following code easier to think about.
    HDirections                  = HDirections + 90;
    HDirections(HDirections>360) = HDirections(HDirections>360) - 360;
    % Split the data by arm occupany, then for each arm check if we are facing CW or CCW %
    % armsDims(n,:) gives the upper and lower limits into arm n, arms start at the TL corner and run clockwise (usual xy->pol conversion is anti-clock, but DACQ coords have inverted y-axis).
    if strcmp(obj.type,'dacq')
        dirDefs                 = [0, 90, 180, -90]; 
    elseif strcmp(obj.type,'npix')
        dirDefs                 = [0, -90, 180, 90];  % These define the amount to be added to the actual dirs, before classifying CW and CCW. One for each arm.
    end
    dirInd                      = zeros(size(HDirections));
    for i = 1:4
        tempDir                 = HDirections;
        ind                     = linPos > armDims(i,1) & linPos <= armDims(i,2);
        tempDir(~ind)           = nan;
        tempDir                 = mod( tempDir+dirDefs(i), 360 );
        dirInd(tempDir<=180)    = 1;  % 1=CW
        dirInd(tempDir>180)     = 2;  % 2=CCW
    end
    
    % We will also define the 'main arm' dwell, to easily exclude corners from the direction definitions data (direction is unreliable, behaviour might be different) %
    % 2017-10-02: this isn't actually being used at the moment, but remains as it seems it might be useful.
    %     cornerLengths = (prms.cornerExtent./100) .* (armDims(:,2)-armDims(:,1));
    %     cornerLimits  = armDims + [cornerLengths, -cornerLengths];
    %     isMainArmPosInd = false(size(linPos));     whichMainArmPosInd = nan(size(linPos));
    %     for i=1:4
    %         ind = linPos>cornerLimits(i,1) & linPos<cornerLimits(i,2);
    %         isMainArmPosInd( ind )    = true;
    %         whichMainArmPosInd( ind ) = i;
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. Remove jumps in light-defined direction data. %
    % The strategy here is to test whether jumps (short epochs of one dir) are going against the overall direction of position
    % movement. So, if we find a jump where the smoothed path gradient is great than a threshold in the opposite direction,
    % remove the jump by reversing the direction.
    % 3a. Smooth the path (to look for larger scale dir trends) %
    linPosRad  = ( linPos ./ max(armDims(:)) ) .* (2*pi);
    linPosUW   = ( unwrap(linPosRad) ./ (2*pi) ) .* max(armDims(:));
    k          = fspecial('gaussian',[ round(obj.mapParams.linpos.filtSigmaForRunDir * obj.trialMetaData(trialInd).posFs * 3), 1], round(obj.mapParams.linpos.filtSigmaForRunDir * obj.trialMetaData(trialInd).posFs) );
    linPosUWSm = imfilter(linPosUW,k,'replicate');
    % 3b. Now remove the jumps %
    stateChangeInd    = diff( dirInd == 1 ) ~= 0;
    stateChangeNumInd = find(stateChangeInd)';
    epochLengths      = [stateChangeNumInd(1),  diff( stateChangeNumInd )];   % List of lengths (in samples)
    isEpochJumpNumInd = find( epochLengths < (obj.mapParams.linpos.durThrJump * obj.trialMetaData(trialInd).posFs) );
    isEpochJumpNumInd = isEpochJumpNumInd( isEpochJumpNumInd~=1 );  % Don't attempt to remove the first epoch, even if it is small.
    for j=isEpochJumpNumInd
        ind = (stateChangeNumInd(j-1)+1) : stateChangeNumInd(j);
        gr  = mean(diff(linPosUWSm(ind)));  % gr is the gradient of the smoothed linear positions in the jump window.
        gr  = (gr * obj.trialMetaData(trialInd).posFs) / (obj.trialMetaData(trialInd).ppm / 100);  % Change units to cm/s
        if dirInd(ind(1)) == 2  && gr > obj.mapParams.linpos.gradThrForJumpSmooth
            dirInd(ind) = 1;
        elseif dirInd(ind(1)) == 1  && gr < -obj.mapParams.linpos.gradThrForJumpSmooth
            dirInd(ind) = 2;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. Define coherent runs (of over a duration threshold) %
    if obj.mapParams.linpos.durThrCohRun > 0
        stateChangeInd    = diff( dirInd == 1 ) ~= 0;
        stateChangeNumInd = find(stateChangeInd)';
        epochLengths      = [stateChangeNumInd(1),  diff( stateChangeNumInd )];   % List of lengths (in samples)
        isCohRunNumInd    = find( epochLengths > ( obj.mapParams.linpos.durThrCohRun * obj.trialMetaData(trialInd).posFs ) );
%             cohRunInd         = zeros(size(xy,1),1);
        cohRunInd         = false(size(xy,1),1);
        
        for j=isCohRunNumInd
            %need to account or case when very first run is coherent
            if j == 1
                ind = 1:stateChangeNumInd(1);
            else
                ind = (stateChangeNumInd(j-1)+1) : stateChangeNumInd(j);
            end           
%             if dirInd(ind(1)) == 1
%                 cohRunInd(ind) = 1;
%             elseif dirInd(ind(1)) == 2
%                 cohRunInd(ind) = 2;
%             end

            if dirInd(ind(1)) ~= 0
                cohRunInd(ind) = true;
            end
        end
        dirInd(~cohRunInd) = 0; % remove non-coherent runs from data
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif contains(lower(obj.trialMetaData(trialInd).trialType),'lin')
    
    if isempty(obj.mapParams.linpos.runDimension)
        [ ~, obj.mapParams.linpos.runDimension ] = max(range(xy)); % estimate run dimension
    end
    
    linPos = xy( :, obj.mapParams.linpos.runDimension );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find ends of track %
    pHist  = histcounts( linPos, 0:1:ceil(max(linPos)) );
    goodP  = pHist >= obj.mapParams.linpos.minDwellForEdge * obj.trialMetaData(trialInd).posFs;
    trEnds = [find(goodP,1,'first'); find(goodP,1,'last')];
    % remove data form beyond ends of track %
    linPos(linPos < trEnds(1))   = NaN;
    linPos(linPos > trEnds(2))   = NaN;
    % Scale %
    linPos                       = linPos - trEnds(1);
    SF                           = trackLength / diff(trEnds);
    linPos                       = linPos .* SF;
    linPos(linPos > trackLength) = trackLength;    % Have checked that when this happens, it is a many decimal places rounding error.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make an index of which direction on the track the rat was running %
    if obj.mapParams.linpos.runDimension == 1
        armDirs = [0 180];   % Remember that the dir data is always in physics convention in matlab (0deg=east).
    else
        armDirs = [90 270];
    end
    dirInd = zeros( size(HDirections) );
    for i = 1:2
        ang2ArmDir    = circ_dist( circ_ang2rad(HDirections), circ_ang2rad(armDirs(i)) );
        ind           = abs(ang2ArmDir) <= obj.mapParams.linpos.dirTolerance;
        dirInd( ind ) = i;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if obj.mapParams.linpos.durThrCohRun > 0
        % Define coherent runs (of over a duration threshold) %
%         cohRunInd         = zeros(size(xy,1),1);
        cohRunInd             = false(size(xy,1),1);
        for i = 1:2
            stateChangeInd    = diff( [0; dirInd == i; 0] );
            stateChangeNumInd = find( stateChangeInd );
            epochLengths      = diff( stateChangeNumInd );       % List of lengths (in samples)
            isCohRunNumInd    = find( epochLengths > (oobj.mapParams.linposptions.durThrCohRun * obj.trialMetaData(trialInd).posFs) );
            for j=1:length(isCohRunNumInd)
                ind = (stateChangeNumInd( isCohRunNumInd(j)) ) : (stateChangeNumInd( isCohRunNumInd(j)+1 ) - 1);
%                 if dirInd(ind(1)) == 1
%                     cohRunInd(ind) = 1;
%                 elseif dirInd(ind(1)) == 2
%                     cohRunInd(ind) = 2;
%                 end
                if dirInd(ind(1)) ~= 0
                    cohRunInd(ind) = true;
                end
            end
        end
        dirInd(~cohRunInd) = 0; % remove non-coherent runs from data
    end
end

end

%% sub function to deal with square track data
function [linPos, armDims] = radialLineariseSquare_v2(XY,envEdges,nPixPerArm)
% Linearise positions from square track data by matching pos samples to 
% bins of an idealised rectangle/square of same size. 
%
% Usage: [linPos,linDir] = radialLineariseSquare_v2(XY,envEdges,nPixPerArm)
%
% Inputs:   posXY (nPosSamples,2) - XY pos coordinates
%           envEdges (2,2)        - edges of environment in x and y (in pixels)
%           nPixPerArm (1)        - how many pix should arms be given ppm and
%                                   actual rectangle size
%
% Outputs:  linPos  - linearised positions
%           armDims - 
%
% TW/LM 2017; LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make the idealised rectangle
% Make an (nPoint,XY) array that descibes an idealised rectangle.
envSize                                 = envEdges(2,:) - envEdges(1,:) + 1;  % +1 so that size is inclusize of both edges, and matches length of vx and vy vectors.
vx                                      = linspace(envEdges(1,1), envEdges(2,1), nPixPerArm)';
vy                                      = linspace(envEdges(1,2), envEdges(2,2), nPixPerArm)';
% Goes CCW round square, starting at TL corner at 1,1.
sqXY( 1:length(vx), 1:2)                = [vx, ones(size(vx)).*envEdges(1,2)];          % south side,
armDims(1,:)                            = [1, length(vx)];
%direction right to left
sqXY( armDims(1,2)+(1:length(vy)), 1:2) = [ones(size(vy)).*envEdges(2,1), vy ];         % West side,
armDims(2,:)                            = [armDims(1,2), armDims(1,2)+length(vy)];
% %direction up to down
sqXY( armDims(2,2)+(1:length(vx)), 1:2) = [flipud(vx), ones(size(vx)).*envEdges(2,2)];  % North side,
armDims(3,:)                            = [armDims(2,2), armDims(2,2)+length(vx)];
% %direction left to right.
sqXY( armDims(3,2)+(1:length(vy)), 1:2) = [ones(size(vy)).*envEdges(1,1), flipud(vy) ]; % East side,
armDims(4,:)                            = [armDims(3,2), armDims(3,2)+length(vy)];

% Centre the square at 0,0 and convert to points in sqXY to radial coords.
sqXY      = bsxfun( @minus, sqXY, (envSize./2)+envEdges(1,:) );
[sqTh,~]  = cart2pol(sqXY(:,1), sqXY(:,2));
sqTh      = mod(sqTh,2*pi);
% Likewise, centre the real data at 0,0 and convert to radial coords.
XY        = bsxfun( @minus, XY, (envSize./2)+envEdges(1,:) );
[posTh,~] = cart2pol(XY(:,1), XY(:,2));
posTh     = mod(posTh,2*pi);

% get radial bin indices
% Need to make sure that the theta coords of the square monotononically
% increase circularly around the square. We need to set theta start (which
% will be in south/west corner) as smallest possible value (in square it 
% would always be 1.25*pi, but in rectangle it's slightly different values)
lowerLimTh   = sqTh(1); % Here we are defining the break point as being the BL corner.
ind          = sqTh<lowerLimTh;
sqTh(ind)    = sqTh(ind) + 2*pi; %bump up coordinates
sqTh         = sqTh(1:end-1); %remove last value as it is == first value 
%also bumb thetas in actual pos vector
ind          = posTh<lowerLimTh;
posTh(ind)   = posTh(ind) + 2*pi;
% Use histc to convert the real data radial positions into an index of
% where they fall with respect to the idealiased square bin limits.
[~,~,linPos] = histcounts(posTh, sqTh);

end



