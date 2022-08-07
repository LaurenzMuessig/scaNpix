function [gridness, Props] = gridprops(autoCorr,varargin)
% Using autocorrelogram calculates grid wavelength, orientation, gridness and x,y
% offset of peaks
% package: scanpix.analysis
%
%       [gridness] = scanpix.analysis.gridprops( autoCorr );
%       [gridness, Props] = scanpix.analysis.gridprops( autoCorr, 'paramName', 'paramValue', .. );
%
% Produces key metrics about a grid cell derived from the spatial autocorrelogram: wavelength,
% orientation, gridness, the mean x,y offset of the six central peaks (if fewer than six peaks are
% found then that number is used to calculate the values). Also returns the coordinates of the peaks
% used to calculate these values [x,y] pairs with origin at centre
%
% Optional input parameters: 
%
%   'peakMode'    - 'point' or 'area'
%   'corrThr'     - Peaks are r > this value
%   'areaThr'     - Peaks are > this many bins (in peakMode=area)
%   'corrThrMode' - 'abs' or 'rel'. Does corrThr refer to absolute  r-value (default) or r relative to central point ('rel', only use for time-win ACs).
%   'radius'      - Radius finding method. 'fieldExtent' (Outermost extent of six peaks), 'est' (mean peak distance *1.25),
%                                          'fix:X' (X=nBin), or 'none' (use whole AC).
%   'fieldExtentMethod'   - 'halfHeight' or 'watershed'. 
%   'cenPkThr     - Threshold used to define extent of central peak (removed from gridness corrs). Relative to central point.
%   'gridMaskSizeThr - If (central peak area)/(outer peak circle area) > centreMaskSizeThr, reject as grid cell. Using this option ..
%                      checks that the AC has a well-defined central peak, otherwise it's not a grid. Default=0.5. Set to [] to inactivate.
%   'crossCorrMode'  - 1 or <0>. Assumes CG is a cross-cell CG. Looks for closest peak to centre, treats this point as the centre of the CG and proceeds
%                      to calculate gridness around this point.  Set 'centreMaskSizeThr' to [], as this code hasn't been fixed for cross-corrs.
%   'verbose'     - Print to screen prms struct
%   'showCorrRing'  - If 1, Display figure with 'doughnut' of autocorr bins that are rotated and correlated to form grid score.
%
%  Fields of additional output properties structure ('Props'):
%
%  .waveLength       NB. Unit for wavelength is bins of autocorr
%  .gridness
%  .orientation
%  .fieldSize
%  .closestPeaksCoord 
%  .maxDistFromCentre 
%  .centralPeakMask
%  .crossCorrClosestPeak   NB. x,y coords, relative to centre pixel. Only valid for cross-corrs.
%
% Note: In the autocorrelogram it is sensible to exclude bins that were constructed with relatively
% small overlap between the ratemap1 and ratemap2 (Hafting excludes bins with an overlap of 20 or
% less). Set these bins to 0 before passing to this function

prms.peakMode = 'point';     % Are peaks local max points ('point'), or contig areas over corrThr ('area')?
prms.corrThr = 0.3;          % For peakMode='point', points must be >corrThr. For peakMode='area', look for contig regions over corrThr.
prms.fieldExtentMethod = 'watershed'; % How do we find the extent of the peaks. 'halfHeight' or 'watershed'.
prms.corrThrMode = 'abs';    % Does corrThr refer to absolute r-value (default) or r relative to central point ('rel', for time-win ACs).
prms.areaThr = 20;           % Min size of contig areas when peakMode='area'. 100(MoserThr) * 1.5^2(MoserBin) / 1.92^2(WillsBin) = 61
prms.closePeakFilter = [0 1 1 1 0; ones(3,5); 0 1 1 1 0];    % When peakMode='point', this filter defines area within which two peaks are counted as one.
prms.radius = 'est';         % How are peaks used to estimate gridness ring. 
prms.cenPkThr=0.5;           % Threshold used to define extent of central peak (removed from gridness corrs). Relative to central point.
prms.centreMaskSizeThr=0.5;    % If (central peak area)/(outer peak circle area) > centreMaskSizeThr, reject as grid cell. Using this option ..
                             %  .. checks that the AC has a well-defined central peak, otherwise it's not a grid. Set to [] to inactivate.
prms.crossCorrMode=0;        % Assumes CG is a cross-cell CG. Looks for closest peak to centre, treats this point as the centre of the CG and proceeds
                             % to calculate gridness. Set 'centreMaskSizeThr' to [], as this code hasn't been fixed for cross-corrs.
prms.verbose=0;              % Print to screen prms struct
prms.showCorrRing=0;         % Display figure with 'doughnut' of autocorr bins that are rotated and correlated to form grid score.


% ---------------------------------------------------------------------------------- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %

if prms.verbose;  disp(prms);   end

if isempty(prms.closePeakFilter); prms.closePeakFilter = [0 1 1 1 0; ones(3,5); 0 1 1 1 0]; end % this is just for compatability with the GUI 

% Preallocate output Props %
gridness=nan;
Props.waveLength=NaN;
Props.gridness=NaN;
Props.orientation=nan(1,3);
Props.fieldSize=NaN;
Props.closestPeaksCoord=NaN; 
Props.maxDistFromCentre=NaN;
Props.centralPeakMask=NaN;
Props.crossCorrClosestPeak=NaN;

% Protect against all(autocorr==nan) (when rate map is 0Hz) %
if all(isnan(autoCorr));   return;     end

% Pre-treat the AC to remove NaNs  %
autoCorr(isnan(autoCorr)) = -1;

% Find Peaks %
if strcmp(prms.peakMode,'point')
    % Find peaks by looking for regional maxima. 
    autoCorrTemp = autoCorr;
    autoCorrTemp(isnan(autoCorrTemp))=-1;   %Sub nans for -1
    if ~isempty(prms.corrThr)
        if strcmp(prms.corrThrMode,'abs')
            autoCorrTemp(autoCorr<=prms.corrThr) = -1;   % Do not allow local max that have r-value below threshold, absolute.
        elseif strcmp(prms.corrThrMode,'rel')
            autoCorrTemp(autoCorr<=(prms.corrThr*max(max(autoCorr)))) = -1;   % Do not allow local max that have r-value below threshold, relative to peak.
        end
    end
    peaksAutoCorr = imregionalmax(autoCorrTemp);    % Find local maxima
    peaksAutoCorr = imdilate(peaksAutoCorr, prms.closePeakFilter); % TW: dilate peaks, so peaks within 4-6 bins are counted as one.
    [lableMask, dummy]=bwlabel(peaksAutoCorr, 8);
elseif strcmp(prms.peakMode,'area')
    % Peaks are areas of contiguous bins larger than areaThr, with corr
    % greater than corrThr (e.g. Sargolini et al, 2006).
    peaksAutoCorr=autoCorr>prms.corrThr;
    [lableMask, dummy]=bwlabel(peaksAutoCorr, 8);
    areas = regionprops(lableMask,'area');
    for ii=1:length(areas)
        if areas(ii).Area < prms.areaThr
            lableMask(lableMask==ii)=0;
        end
    end
    [lableMask, dummy]=bwlabel(lableMask, 8);    % Re-find surviving peaks
end
% If no peaks not found, reject as grid cell %
if length(unique(lableMask))<=2;   return;    end

    
%In case adjacent points share maxima find the centroid of them - NB returns structure array
%stats(n).Centroid containing for each peak the x,y position of the centroid but y is
%counting down from origin at top left
stats=regionprops(lableMask, 'Centroid');
%NL. [n x 2] pairs of x,y coord for max points
xyCoordMaxBin=round(reshape([stats.Centroid], 2,[])'); %Still x,y pair

% Convert to a new reference frame which as the origin at the centre of the autocorr
% NB autocorr will always have sides with odd number of bins
centralPoint=ceil(size(autoCorr)/2); %m,n pair
xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPoint), [length(xyCoordMaxBin), 1]);

%Calculate distance of peaks from centre point and find seven closest (one will be central peak
%disregard this)
distFromCentre=sum(xyCoordMaxBinCentral.^2,2).^0.5;
[dummy, orderOfClose]= sort(distFromCentre);

%Get id of closest peaks and central peak
centralPeak=orderOfClose(1);
if length(orderOfClose)>=7; closestPeaks=orderOfClose(2:7); %Might be fewer than 7 peaks
else closestPeaks=orderOfClose(2:end);
end

%%%% Cross-correlogram Mode. To measure gridness in a cross-correlogram, find closest peak to centre        %%%%%
%%%% and subtract all other peak positions from this, ie treat this point as the centre of the correlogram. %%%%%
if prms.crossCorrMode
    closestPeakCoordsCentre = xyCoordMaxBinCentral(centralPeak,:);
    Props.crossCorrClosestPeak = closestPeakCoordsCentre;   % Assign this peak to output structure.
    Props.corrCorrOriginalClosestPeaksCoord = xyCoordMaxBin([centralPeak; closestPeaks],:);  % This is a record of where the peaks were before shifting.

    % STEP 1: Offset the CG so that closest peak is at centre
    padAC = nan(size(autoCorr).*3);
    shiftIndexRow = (size(autoCorr,1)+1 : size(autoCorr,1)*2) - closestPeakCoordsCentre(2);
    shiftIndexCol = (size(autoCorr,2)+1 : size(autoCorr,2)*2) - closestPeakCoordsCentre(1);
    padAC(shiftIndexRow,shiftIndexCol) = autoCorr;
    
    % STEP 2: Also subtract closest peak co-ords from all, so that closest becomes centre, all other shift in unison with CG
    xyCoordMaxBin =  xyCoordMaxBin - repmat(closestPeakCoordsCentre, size(xyCoordMaxBin,1), 1) + repmat( fliplr(size(autoCorr)), size(xyCoordMaxBin,1), 1 );
     
    % STEP 3: Re-find closest peaks, same code as previously should work.
    autoCorr=padAC;
    centralPoint=ceil(size(autoCorr)/2); % r,c pair
    xyCoordMaxBinCentral=xyCoordMaxBin-repmat(fliplr(centralPoint), [length(xyCoordMaxBin), 1]);
    %Calculate distance of peaks from centre point and find seven closest (one will be central peak disregard this)
    distFromCentre=sum(xyCoordMaxBinCentral.^2,2).^0.5;
    [dummy, orderOfClose]= sort(distFromCentre);
    % Get id of closest peaks and central peak centralPeak=orderOfClose(1) ;
    centralPeak=orderOfClose(1);
    if length(orderOfClose)>=7; closestPeaks=orderOfClose(2:7); %Might be fewer than 7 peaks
    else closestPeaks=orderOfClose(2:end);
    end
    
end

%x,y pairs in cartesian coords with y counting down from top
Props.closestPeaksCoord=xyCoordMaxBin(closestPeaks,:);

% -------------------------------------------------------------------------------------------------
% --- MEAN X & y OFFSET OF 6 CENTRAL PEAKS --------------------------------------------------------
% -------------------------------------------------------------------------------------------------

meanXyOffset=mean((xyCoordMaxBinCentral(closestPeaks,:).^2).^0.5,1);


% -------------------------------------------------------------------------------------------------
% --- FIND FIELD AROUND EACH PEAK -----------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%As we've just defined peak of each field (when prms.peakMode='point') need to find the extent of the field around it. Define
%this as the area enclosed within the half-height of the peak. Do this by looping through eack peak:
%find areas of SAC above the half height, then find the patch that includes the peak of interest.
% TW: original CB code introduced 'perimeterMask' variable, to help recover peak extent when peaks
%     overlapped. Note that this isn't used in current code (as only use for peak extents is to look 
%     for the furthest from the centre), but the variable is retained.

% Find extent of 6-closest peaks %
if strcmp(prms.peakMode,'point')
    if strcmp(prms.fieldExtentMethod,'halfHeight')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find the half-height around each peak %
        % Preassign output: Two 3d mats [sizeAutoCorr1, sizeAutoCorr2, numberPeaks]
        peakMasks=zeros([size(autoCorr), size(Props.closestPeaksCoord,1)]);
        perimeterMasks=zeros([size(autoCorr), size(Props.closestPeaksCoord,1)]);
        for n=1:size(Props.closestPeaksCoord,1)
            %NB. Have to flip dimensions as these are x,y pairs and need to be m,n
            [peakMasks(:,:,n+1), perimeterMasks(:,:,n+1)]=findPeakExtent(autoCorr, closestPeaks(n), [Props.closestPeaksCoord(n,2), Props.closestPeaksCoord(n,1)]);
        end
        lableMask=max(peakMasks, [], 3);
        perimeterMask=max(perimeterMasks,[],3);
        
    elseif strcmp(prms.fieldExtentMethod,'watershed')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %2a. Find the inverse-drainage-basin for each peak.
        fieldsLabel = watershed(-autoCorr) + 1; %plus 1, because we want to be able to use label as index

        %2b. Work out what threshold to use in each drainage-basin
        nZones = max(fieldsLabel(:));
        closePeakInds = sub2ind(size(autoCorr),Props.closestPeaksCoord(:,2),Props.closestPeaksCoord(:,1));
        closesPeaksNewId = fieldsLabel(closePeakInds);
        thresholds = Inf(nZones,1);
        thresholds(closesPeaksNewId) = autoCorr(closePeakInds)/2;

        %2c. Apply thresholds to get a mask and updated labels
        fieldsMask = autoCorr > thresholds(fieldsLabel);
        fieldsLabel(~fieldsMask) = 0; %note label numbers are the ones from watershed not from finding peaks

    end
end
% Modified, TW. Treat central peak separately. Has a dedicated threshold (prms.cenPkThr). Always look for extent, and keep this in a dedicated variable %
if ~isempty(prms.cenPkThr)
    Props.centralPeakMask=findPeakExtent(autoCorr,centralPeak,centralPoint,prms.cenPkThr);
    Props.centralPeakMask(Props.centralPeakMask>0)=1;
    Props.centralPeakMask=logical(Props.centralPeakMask);
else
    Props.centralPeakMask=false(size(autoCorr));
end

% -------------------------------------------------------------------------------------------------
% --- WAVELENGTH ----------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------

%Calculate wavelength measured in bins
% waveLength=median(distFromCentre(closestPeaks));
Props.waveLength=mean(distFromCentre(closestPeaks));

% -------------------------------------------------------------------------------------------------
% --- ORIENTATION ---------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%Calculate orientation of grid - note orientation is measured counterclockwise from horizontal so
%basically polar coordinates
[th, dummy]=cart2pol(xyCoordMaxBinCentral(closestPeaks,1), -xyCoordMaxBinCentral(closestPeaks,2)); % TW: Modified to take account of inverted y-axis.
th=th(th>=0); %Remove negative values - can do this as peaks are 180deg radially symetrical
if ~isempty(th)
%     Props.orientation=(min(th)/(2*pi))*360;
    tmp = (sort(th)./(2*pi))*360;
    if length(tmp) > 3
        Props.orientation=tmp(1:3)';
    else
        Props.orientation(1:length(th))=tmp;
    end
end
clear r th

% --------------------------------------------------------------------------------------------------
% ---- GRIDNESS ------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
[meshX, meshY]=meshgrid(1:size(autoCorr,2), 1:size(autoCorr,1));
%Change coordinates so that origin is at centre and to facilitate pythag
meshXCentre=(meshX-centralPoint(2)).^2;
meshYCentre=(meshY-centralPoint(1)).^2;
distToCentre = (meshXCentre+meshYCentre).^0.5;

% Find outer radius of the gridness calculation zone (maxDistFromCentre) %
if strcmp(prms.radius,'fieldExtent')
    
    % Method 1: Find extent of the 6 closest peaks and take the max %
    if strcmp(prms.fieldExtentMethod,'halfHeight')
        maxDistTemp=zeros(length(closestPeaks),1);
        for ii=1:length(closestPeaks)
            isPresent = lableMask==closestPeaks(ii);
            %Deal with situations where the mask of one peak has bloted out another entirely
            if ~sum(isPresent(:)), maxDistTemp(ii)=0; continue, end
            maxDistTemp(ii)=max((meshXCentre(lableMask==closestPeaks(ii)) + meshYCentre(lableMask==closestPeaks(ii))).^0.5);
        end
    elseif strcmp(prms.fieldExtentMethod,'watershed')
        maxDistTemp = accumarray(fieldsLabel(fieldsMask),distToCentre(fieldsMask),[],@max);
    end
    Props.maxDistFromCentre=max(maxDistTemp);
    
elseif strcmp(prms.radius,'est') 
    % Method 2: wavelength plus assumed half-height %
    Props.maxDistFromCentre = Props.waveLength * 1.25;
elseif strcmp(prms.radius(1:3),'fix')
    Props.maxDistFromCentre = str2double(prms.radius(5:end));
elseif strcmp(prms.radius,'none')
    Props.maxDistFromCentre = min(floor(size(autoCorr)/2));   % Find max radius fitting in AC
end

%MaxDistFromCentre can not be greater than the distance from the centre to the nearest edge of
%autoCorr - if it is larger reset it to this lower value
if Props.maxDistFromCentre>min(floor(size(autoCorr)/2))
    Props.maxDistFromCentre=min(floor(size(autoCorr)/2));
end

%Create a circular mask using the meshgrids I've created above
gridnessMaskOuter=distToCentre<=Props.maxDistFromCentre; % Mask for points within bound of outer peak
gridnessMask=gridnessMaskOuter & ~(Props.centralPeakMask);                   % Gridness 'ring' is this minus central peak.
if ~isempty(prms.centreMaskSizeThr)
    % Following line says that if the central peak is > prms.centreMaskSizeThr (e.g. 0.5) than the area within
    % the outer peak, reject as possible grid cell. Protects against ACs that
    % have large central peak with a few wispy bits on the edge that are hexagonal.
    % Has additional effect that time-win ACs without well-defined central peak will be
    % rejected, as central peak will be large. (prms.cenPkThr is always relative to max r of peak).
    if sum(sum(Props.centralPeakMask))/sum(sum(gridnessMaskOuter)) > prms.centreMaskSizeThr;   return;   end
end
clear meshX meshY

%%%%%%%%%%%%%% Calculate gridness correlations %%%%%%%%%%%%
selA=autoCorr;
selA(~gridnessMask)=NaN;
if prms.showCorrRing;	figure;  imagesc(selA);  end  % For easy testing of what gridness calcualtion is working with.
rotAmnt = [60, 120, 30, 90, 150];
rotatedCorr = nan(1,5);
for ii=1:length(rotAmnt)
    selB=imrotate(autoCorr,rotAmnt(ii), 'bilinear', 'crop');
    selB(~gridnessMask)=NaN;
    visMask = intersect( find(~isnan(selA)), find(~isnan(selB)) );
    visA=selA(visMask);   visB=selB(visMask);
    rotatedCorr(ii)=corr2(visA,visB);
end

%Finally gridness is the difference between the lowest correlation at 60degand 120deg and the highest correlation at 30, 90 and 150deg
gridness= min(rotatedCorr([1, 2])) - max(rotatedCorr([3, 4, 5])); %This one for quick verison
Props.gridness = gridness;

% -------------------------------------------------------------------------
% --- FIELD SIZE (Hafting 2005)--------------------------------------------
% -------------------------------------------------------------------------
% Size of central peak, using r=0.4 as threshold.
aboveThrMask = bwlabel((autoCorr>0.4),8); 
aboveThrMask( aboveThrMask ~= aboveThrMask(centralPoint(1),centralPoint(2)) ) = 0;
aboveThrMask( aboveThrMask == aboveThrMask(centralPoint(1),centralPoint(2)) ) = 1;
aboveThrMask=logical(aboveThrMask);
% If central peak bleeds to outside of AC, don't count field size %
unVisMask=zeros(size(autoCorr));
unVisMask(isnan(autoCorr))=1;
tempLabel=bwlabel( unVisMask|aboveThrMask ,8);
if length(unique(tempLabel))~=2
    Props.fieldSize=sum(sum(aboveThrMask)); 
end


% -------------------------------------------------------------------------------------------------
% --- INLINE FUNCTIONS ----------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------

function [peakMask, perimeterMask]=findPeakExtent(autoCorr,peakId,peakCoord,varargin)
%Finds extent of field that belongs to each peak - defined as area in half-height and also
%perimieter. NB. peakCoord must by m,n pair in normal matrix coords
% The last argument defines the threshold level. This can be omitted, and is 0.5 by default.
if isempty(varargin);   pkThr=0.5;  else   pkThr=varargin{1};   end
peakMask=zeros(size(autoCorr));
perimeterMask=zeros(size(autoCorr));
aboveHalfHeightMask=bwlabel( autoCorr > (autoCorr(peakCoord(1),peakCoord(2))*pkThr ),8); 
peakIdTemp=aboveHalfHeightMask(peakCoord(1),peakCoord(2));
peakMask(aboveHalfHeightMask==peakIdTemp)=peakId;
perimeterMask(bwperim(aboveHalfHeightMask==peakIdTemp))=peakId;






% % ------------------------------------------------------------------------------------------------
% % --- CB MEASURE OF GRID REGULARITY --------------------------------------------------------------
% % ------------------------------------------------------------------------------------------------
% %Basically get the distance to the three peaks of the SAC - find difference between the X peak (one
% %closest to the X axis) and the Y peak (one closest to Y axis)
% %Have changes this now to be dist to X peak divided by dist to Y peak. Gives a value that is
% %comparable across cells with different wavelengths
% [th, r]=cart2pol(xyCoordMaxBinCentral(closestPeaks), xyCoordMaxBinCentral((closestPeaks)+length(xyCoordMaxBinCentral)));
% [sortedTh,sortInd]=sort(th);
% sortedR=r(sortInd);
% peakClosestToX=find(abs(sortedTh)==min(abs(sortedTh)),1);
% peakClosestToY=find(abs(sortedTh-(pi/2))==min(abs(sortedTh-(pi/2))),1);
% distToSACXPeak=sortedR(peakClosestToX);
% distToSACYPeak=sortedR(peakClosestToY);
% % distToSACXYPeak=[distToSACYPeak, distToSACXPeak];
% distToSACXYPeak=distToSACXPeak/distToSACYPeak;
% 
% %
% cbGridness=(sortedR(peakClosestToX)-sortedR(peakClosestToY))/waveLength;


% -------------------------------------------------------------------------
% --- CHECK HEXAGONALITY OF PEAKS -----------------------------------------
% -------------------------------------------------------------------------
% TW: not in use. Routine not good enough at finding true peaks - rejects
% too many grid cells.
% Reject grid if peaks not hexagonal enough. Not used by default. Caller
% must enter 'hexFilt' input arg.
% if ~isempty(prms.hexFilt)
%     degDiff = diff( sort( (th/(2*pi))*360 ) );
%     if any( degDiff<(60-prms.hexFilt) | degDiff>(60+prms.hexFilt) )
%     	waveLength = nan;   gridness = nan;   orientation = nan;
%         fieldSize = nan;   closestPeaksCoord = nan;   maxDistFromCentre = nan;
%         return
%     end
% end
