function loadPosNPix(obj, trialIterator)
% loadPos - load position data for neuropixel data
%
% Syntax:  loadPos(obj, trialIterator)
%
% Inputs:
%    obj           - ephys class object ('npix')
%    trialIterator - numeric index for trial to be loaded
%
% Outputs:
%
% See also: 
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialIterator (1,1) {mustBeNumeric}
end

%%
fprintf('Loading pos data for %s .......... ', obj.trialNames{trialIterator});

%% process data

% open pos data the format is [frame count, greenXY, redXY winSzX, winSzY, timeStamp possibly other Data ]
fName = dir(fullfile(obj.dataPath{trialIterator},'trackingData', '*.csv'));
if isempty(fName)
    warning(['scaNpix::npixUtils::loadPosNPix:Can''t find csv file in ' obj.dataPath{trialIterator} '. Come on mate.']);
    return;
end

fID = fopen(fullfile(fName.folder,fName.name),'rt');
header = textscan(fID,'%s',1);
nColumns = length(strsplit(header{1}{1},','));
fmt = '%u%f%f%f%f%u%u%f';
% allow for any n of additonal fields from Bonsai output
if nColumns > 8; fmt = [fmt repmat('%u',nColumns-8,1)]; end

csvData = textscan(fID,fmt,'HeaderLines',1,'delimiter',',');
fclose(fID);

% led data 
if strcmp(obj.trialMetaData(trialIterator).LEDfront,'green')
    led          = [csvData{2}, csvData{3}]; % xy coords
    led(:,:,2)   = [csvData{4}, csvData{5}]; % xy coords
else
    led          = [csvData{4}, csvData{5}]; % xy coords
    led(:,:,2)   = [csvData{2}, csvData{3}]; % xy coords
end
led(led==0)      = NaN;

% sample Times
timeStamps       = csvData{8};
sampleT          = scanpix.npixUtils.convertPointGreyCamTimeStamps(timeStamps); % starts @ 0

% in case logging point grey data was corrupt
if all(sampleT == 0)
    sampleT    = (0:length(led)-1)' * 1/obj.trialMetaData(trialIterator).posFs; % pretend we have perfect sampling
    frameCount = [1;(length(led):-1:2)'+10e2]; % make a mock frame count that is corrupt from sample 1 onwards so we can use the fix in 'fixFrameCount' (in-line func.)
    obj.trialMetaData(trialIterator).BonsaiCorruptFlag = true;
    warning('scaNpix::npixUtils::loadPosNPix:Point Grey data corrupt!');
else
    frameCount = csvData{1} - csvData{1}(1) + 1;
    obj.trialMetaData(trialIterator).BonsaiCorruptFlag = false;
end

% deal with problems between data streams (inline function)
[frameCount, sampleT] = fixFrameCounts(obj,trialIterator,frameCount,sampleT);

% deal with missing frames (if any) - this currently doesn't take into account if 1st frame(s) would be missing, but I am not sure this would
% actually ever happen (as 1st frame should always be triggered fine)
% first check if there are any...
missFrames       = find(~ismember(1:frameCount(end),frameCount));
nMissFrames      = length(missFrames);
if ~isempty(missFrames)
    fprintf('Note: There are %i missing frames in tracking data for %s.\n', nMissFrames, obj.trialMetaData(trialIterator).filename);
    
    temp                     = zeros(length(led)+nMissFrames, 2, obj.trialMetaData(trialIterator).nLEDs);
    temp(missFrames,:,:)     = nan;
    temp(temp(:,1)==0,:,:)   = led;
    led                      = temp;
    
    % interpolate sample times
    interp_sampleT           = interp1(double(frameCount), sampleT, missFrames);
    temp2                    = zeros(length(led),1);
    temp2(missFrames,1)      = interp_sampleT;
    temp2(temp2(:,1) == 0,1) = sampleT;
    sampleT                  = temp2;   
    %
    obj.trialMetaData(trialIterator).log.missingFramesPosStream = nMissFrames;
end

ppm = nan(2,1);
if isempty(regexp(obj.trialMetaData(trialIterator).trialType,'circle','once')) && size(obj.trialMetaData(trialIterator).envBorderCoords,2) ~= 3; circleFlag = false; else; circleFlag = true; end
% estimate ppm
if isempty(obj.trialMetaData(trialIterator).envBorderCoords)
    envSzPix  = [double(csvData{6}(1)) double(csvData{7}(1))];
    ppm(:)    = mean(envSzPix ./ (obj.trialMetaData(trialIterator).envSize ./ 100) );
else
    % this case should be default
    if ~circleFlag
        % recover all corner coords from 2 points - this should be independent of box misalignment with cam window
        knownDist = sqrt( (obj.trialMetaData(trialIterator).envBorderCoords(1,1)-obj.trialMetaData(trialIterator).envBorderCoords(1,2))^2 + (obj.trialMetaData(trialIterator).envBorderCoords(2,1)-obj.trialMetaData(trialIterator).envBorderCoords(2,2))^2 );
        ppm(:) = round( mean( knownDist ./ (sqrt(sum(obj.trialMetaData(trialIterator).envSize.^2)) ./ 100) ) );
        % full set
        obj.trialMetaData(trialIterator).envBorderCoords = scanpix.helpers.findBoxCorners(obj.trialMetaData(trialIterator).envBorderCoords(:,1),ppm(1)*(obj.trialMetaData(trialIterator).envSize(1)/100), obj.trialMetaData(trialIterator).envBorderCoords(:,2),ppm(1)*(obj.trialMetaData(trialIterator).envSize(2)/100));
        % % now align env coords with the camera window
        % for i = 1:2
        %     led(:,:,i) = scanpix.helpers.rotatePoints(led(:,:,i),[obj.trialMetaData(trialIterator).envBorderCoords(1,1),obj.trialMetaData(trialIterator).envBorderCoords(1,2);obj.trialMetaData(trialIterator).envBorderCoords(2,1),obj.trialMetaData(trialIterator).envBorderCoords(2,2)]);
        % end
        %                     envSzPix  = [abs(obj.trialMetaData(trialIterator).envBorderCoords(1,1)-obj.trialMetaData(trialIterator).envBorderCoords(1,2)), abs(obj.trialMetaData(trialIterator).envBorderCoords(1,3)-obj.trialMetaData(trialIterator).envBorderCoords(2,3))];
    else
        [xCenter, yCenter, radius, ~] = scanpix.fxchange.circlefit(obj.trialMetaData(trialIterator).envBorderCoords(1,:), obj.trialMetaData(trialIterator).envBorderCoords(2,:));
        envSzPix = [2*radius 2*radius];
        ppm(:) = round( mean( envSzPix ./ (obj.trialMetaData(trialIterator).envSize ./ 100) ) );
    end
    %                 ppm(:) = round( mean( envSzPix ./ (obj.trialMetaData(trialIterator).envSize ./ 100) ) );
end

%% post process basically as scanpix.dacqUtils.postprocess_data_v2
% scale data to standard ppm if desired
if ~isempty(obj.params('ScalePos2PPM'))
    scaleFact = (obj.params('ScalePos2PPM')/ppm(1));
    led = floor(led .* scaleFact);
    ppm(1) = obj.params('ScalePos2PPM');
    % obj.trialMetaData(trialIterator).objectPos = obj.trialMetaData(trialIterator).objectPos .* scaleFact;
    obj.trialMetaData(trialIterator).envBorderCoords = obj.trialMetaData(trialIterator).envBorderCoords .* scaleFact;
    if circleFlag
        [xCenter, yCenter, radius] = deal(xCenter*scaleFact,yCenter*scaleFact,radius*scaleFact);
    end
    obj.trialMetaData(trialIterator).PosIsScaled = true;
else
    obj.trialMetaData(trialIterator).PosIsScaled = false;
end

% remove tracking errors that fall outside box
for i = 1:2
    % env borders
    borderTolerancePix = ppm(1)/100*2.5; % we'll assume 1 standard rate map bin tolerance
    if ~circleFlag
        envSzInd = led(:,1,i) < min(obj.trialMetaData(trialIterator).envBorderCoords(1,:))-borderTolerancePix | led(:,1,i) > max(obj.trialMetaData(trialIterator).envBorderCoords(1,:))+borderTolerancePix | led(:,2,i) < min(obj.trialMetaData(trialIterator).envBorderCoords(2,:))-borderTolerancePix | led(:,2,i) > max(obj.trialMetaData(trialIterator).envBorderCoords(2,:))+borderTolerancePix;
    else
        envSzInd = (led(:,1,i) - xCenter).^2 + (led(:,2,i) - yCenter).^2 > (radius+borderTolerancePix)^2; % points outside of environment
    end
    % filter out 
    led(envSzInd,:,i) = NaN;
end

% fix positions (inline subfunction)
% use median frame interval - camera time stamps can be corrupt (see fixFrameCounts) and a mean would be dominated by the corrupt jumps
led = fixPositions(led, median(diff(sampleT)), ppm(1), obj, trialIterator );

% smooth
kernel = ones( ceil(obj.params('posSmooth') * obj.params('posFs')), 1)./ ceil( obj.params('posSmooth') * obj.params('posFs') ); % as per Ephys standard - 400ms boxcar filter
% Smooth lights individually, then get direction.
% doing the smoothing with convolution directly rather than imfilter will prevent spreading of NaNs in the data
smLight = nan(size(led));
for i = 1:2
    smLight(:,:,i) = scanpix.helpers.smoothWithNaNs(led(:,:,i),kernel);
end

% some sanity checks for the data loading
scanpix.npixUtils.dataLoadingReport(length(sampleT),length(obj.spikeData.sampleT{trialIterator}),obj.trialMetaData(trialIterator).BonsaiCorruptFlag);
obj.trialMetaData(trialIterator).log.SyncMismatchPosAP = length(sampleT)-length(obj.spikeData.sampleT{trialIterator});

% align pos data with sync data
endIdxNPix                                = min( [ length(obj.spikeData.sampleT{trialIterator}), find(obj.spikeData.sampleT{trialIterator} < obj.trialMetaData(trialIterator).duration,1,'last') + 1]);
obj.spikeData.sampleT{trialIterator}      = obj.spikeData.sampleT{trialIterator}(1:endIdxNPix);
smLight                                   = smLight(1:endIdxNPix,:,:);

% interpolate positions to pos fs exactly - this will speed up map making significantly 
if obj.trialMetaData(trialIterator).log.InterpPos2PosFs 
    sampleTimes = obj.spikeData.sampleT{trialIterator};
    % newT = (0:1/obj.params('posFs'):(length(sampleTimes)-1)*(1/obj.params('posFs')))'; %
    newT        = linspace(0,sampleTimes(end),length(sampleTimes))'; %
    realFs      = length(sampleTimes) / sampleTimes(end);

    if size(smLight,1) - length(sampleTimes) == 1
        newTint            = newT(2) - newT(1);
        % newT(end+1) = newT(end) + 1/obj.params('posFs');
        % sampleTimes(end+1) = sampleTimes(end) + 1/obj.params('posFs');
        newT(end+1)        = newT(end) + newTint;
        sampleTimes(end+1) = sampleTimes(end) + newTint;
    elseif size(smLight,1) - length(sampleTimes) > 1
        error('scaNpix::npixUtils::loadPosNPix:Something went wrong here. Mismatch of n of pos frames and sync pulses for %s!',obj.trialnames{trialIterator});
    end
    obj.trialMetaData(trialIterator).log.InterpPosSampleTimes = newT;
    obj.trialMetaData(trialIterator).log.InterpPosFs          = realFs;
    obj.trialMetaData(trialIterator).posFs                    = realFs;
    %
    for i = 1:size(smLight,3)
        for j = 1:2
            smLight(:,j,i) = interp1(sampleTimes, smLight(:,j,i), newT);
        end
    end
end

% Get position from smoothed individual lights %% 
wghtLightFront = (1-obj.params('posHead'));
wghtLightBack  = obj.params('posHead');
xy             = smLight(:,:,1) .* wghtLightFront + smLight(:,:,2) .* wghtLightBack;  %

% single LED fallback - where only one LED was tracked (other one lost for > 1s, see fixPositions), use that LED rather
% than discarding the position (head direction stays NaN there, as it needs both LEDs)
[xy, singleLEDInd] = singleLEDFallback(xy, smLight, obj.trialMetaData(trialIterator).posFs);
% where both LEDs were lost, interpolate the position across gaps of up to 'maxPosInterpolate' (s); longer gaps stay NaN
[xy, posInterpInd] = interpShortGaps(xy, round(obj.params('maxPosInterpolate') * obj.trialMetaData(trialIterator).posFs));

obj.trialMetaData(trialIterator).log.posSingleLEDInd  = singleLEDInd;            % samples where position is from one LED only
obj.trialMetaData(trialIterator).log.posSingleLEDFrac = mean(singleLEDInd);
obj.trialMetaData(trialIterator).log.posInterpInd     = posInterpInd;            % samples where position is interpolated (both LEDs lost)
obj.trialMetaData(trialIterator).log.posInterpFrac    = mean(posInterpInd);
obj.trialMetaData(trialIterator).log.posValidFrac     = mean(~isnan(xy(:,1)));   % final fraction of valid positions
if any(singleLEDInd) || any(posInterpInd) || any(isnan(xy(:,1)))
    fprintf('\nNote: positions from a single LED: %.1f%%, interpolated (both LEDs lost): %.1f%%, NaN (both LEDs lost > %.1fs): %.1f%%.\n', 100*mean(singleLEDInd), 100*mean(posInterpInd), obj.params('maxPosInterpolate'), 100*mean(isnan(xy(:,1))));
end

% get direction data
correction     = obj.trialMetaData(trialIterator).LEDorientation(1); %To correct for light pos relative to rat subtract angle of large light
dirData        = mod((180/pi) .* atan2(smLight(:,2,1)-smLight(:,2,2), smLight(:,1,1)-smLight(:,1,2)) - correction, 360); % NaN unless both LEDs tracked (or LED gap <= 1s)
obj.trialMetaData(trialIterator).log.dirValidFrac = mean(~isnan(dirData));

% pos data output
obj.posData(1).XYraw{trialIterator}       = led;
obj.posData(1).XY{trialIterator}          = xy; %[floor(xy(:,1)) + 1, floor(xy(:,2)) + 1];
obj.posData(1).sampleT{trialIterator}     = sampleT(1:endIdxNPix); % this is redundant as we don't want to use the sample times from the PG camera
obj.posData(1).direction{trialIterator}   = dirData;

obj.trialMetaData(trialIterator).ppm      = ppm(1);
obj.trialMetaData(trialIterator).ppm_org  = ppm(2);

% scale position to physical size of environment (only if requested with 'scalePos2Env'; params from older versions lack the key -> off)
if isKey(obj.params,'scalePos2Env') && obj.params('scalePos2Env') && ~obj.params('scalePos2CamWin') && ~isempty(obj.trialMetaData(trialIterator).envSize )
    boxExt = obj.trialMetaData(trialIterator).envSize / 100 * obj.trialMetaData(trialIterator).ppm;
    scanpix.maps.scalePosition(obj, trialIterator,'envSzPix', boxExt,'circleFlag',circleFlag);
else
    obj.trialMetaData(trialIterator).PosIsFitToEnv = {false,[]};
end

% running speed
pathDists                                  = sqrt( diff(xy(:,1)).^2 + diff(xy(:,2)).^2 ) ./ ppm(1) .* 100; % distances in cm
if obj.trialMetaData(trialIterator).log.InterpPos2PosFs
    obj.posData(1).speed{trialIterator}    = pathDists ./ (1/realFs); % cm/s
else
    obj.posData(1).speed{trialIterator}    = pathDists ./ diff(obj.spikeData(1).sampleT{trialIterator}); % cm/s
end
obj.posData(1).speed{trialIterator}(end+1) = obj.posData(1).speed{trialIterator}(end);

fprintf('  DONE!\n');

end

%%%%%%%%%%%%%%%%%%%%%% INLINE FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xy, singleLEDInd] = singleLEDFallback(xy, smLight, posFs)
% fill positions where only one LED is valid. The combined position differs from a single LED by an offset (fraction of
% the LED separation vector, e.g. half of it for posHead = 0.5) that rotates with the head, so is unknown inside the gap.
% We take the offset measured at the edges of each gap (last/first sample with both LEDs) and let it fade out with time
% from the edge (exp. decay, tau = 0.5s). This keeps the position continuous at the gap edges (no speed artefacts) without
% guessing the head direction for long: further inside the gap the single LED is used as is, i.e. the error is ~ the
% offset itself. (Linearly interpolating the edge offsets across the gap was tested too - better on average for short
% gaps, but larger errors (up to ~2x offset) whenever the head turned during the gap.)

tau          = 0.5; % s
nSamp        = size(xy,1);
ledOK        = squeeze(~isnan(smLight(:,1,:))); % nSamp x 2
bothOK       = all(ledOK,2);
singleLEDInd = false(nSamp,1);

for k = 1:2 % k = LED that is still tracked
    onlyK = ledOK(:,k) & ~ledOK(:,3-k);
    if ~any(onlyK); continue; end
    offset   = xy - smLight(:,:,k); % valid where both LEDs are tracked
    d        = diff([false; onlyK; false]);
    runStart = find(d == 1);
    runEnd   = find(d == -1) - 1;
    for r = 1:length(runStart)
        ind = (runStart(r):runEnd(r))';
        pre = runStart(r) - 1;
        pst = runEnd(r) + 1;
        hasPre = pre >= 1 && bothOK(pre);
        hasPst = pst <= nSamp && bothOK(pst);
        % weight of edge offsets, decaying with time from the respective edge
        wPre   = hasPre .* exp( -((ind - pre) ./ posFs) ./ tau );
        wPst   = hasPst .* exp( -((pst - ind) ./ posFs) ./ tau );
        wSum   = max(1, wPre + wPst); % normalise in short gaps, where both edges contribute fully
        offsetInGap = zeros(length(ind), 2);
        if hasPre; offsetInGap = offsetInGap + (wPre ./ wSum) .* offset(pre,:); end
        if hasPst; offsetInGap = offsetInGap + (wPst ./ wSum) .* offset(pst,:); end
        xy(ind,:) = smLight(ind,:,k) + offsetInGap;
    end
    singleLEDInd = singleLEDInd | onlyK;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ledPos = fixPositions(ledPos,frameInt,ppm,obj,trialIterator)
%%
% remove probable tracking errors (i.e. too fast) as per usual
for i = 1:2 
    ok_pos   = find( ~isnan(ledPos(:,1,i)) );
    prev_pos = ok_pos(1);
    for j = 2:length(ok_pos)
        % Get speed of shift from prev_pos
        currSpeed = (sqrt((ledPos(ok_pos(j),1,i)-ledPos(prev_pos,1,i))^2+(ledPos(ok_pos(j),2,i)-ledPos(prev_pos,2,i))^2) / ppm)/ ((ok_pos(j)-prev_pos) * frameInt) ;
        if currSpeed > obj.params('posMaxSpeed')
            ledPos(ok_pos(j),:,i) = NaN;
        else
            prev_pos              = ok_pos(j);
        end
    end
end

%%
% LED pair consistency - where both LEDs are tracked, their separation can't be much larger than the distance between the
% LEDs on the headstage, so a separation > 2x the median separation means one LED was mistracked (e.g. a reflection). We
% remove the LED that jumped, i.e. the one further away from its own (median) position in the surrounding frames. Samples
% where only one LED is tracked are left alone
nSamp   = size(ledPos,1);
bothOK  = all(~isnan(squeeze(ledPos(:,1,:))),2);
LEDsep  = sqrt( sum( (ledPos(:,:,1) - ledPos(:,:,2)).^2, 2) );
badSep  = bothOK & LEDsep > 2 * median(LEDsep(bothOK));
[rem1, rem2] = deal(false(nSamp,1));
if any(badSep)
    dev = nan(nSamp,2);
    for i = 1:2
        dev(:,i)  = sqrt( sum( (ledPos(:,:,i) - movmedian(ledPos(:,:,i), 11, 1, 'omitnan')).^2, 2) ); % deviation from own +/-5 frame median
        % an LED with few valid neighbours (e.g. appearing right after a gap) can't be judged by its own median and is the
        % suspect one (isolated detections are more likely errors)
        nNeighb   = movsum(double(~isnan(ledPos(:,1,i))), 11) - double(~isnan(ledPos(:,1,i)));
        dev(nNeighb < 4, i) = Inf;
    end
    [~, worseLED] = max([sum(isnan(ledPos(:,1,1))), sum(isnan(ledPos(:,1,2)))]); % tie breaker: LED that is tracked worse overall
    rem1 = badSep & (dev(:,1) > dev(:,2) | (dev(:,1) == dev(:,2) & worseLED == 1));
    rem2 = badSep & ~rem1;
    ledPos(rem1,:,1) = NaN;
    ledPos(rem2,:,2) = NaN;
end
obj.trialMetaData(trialIterator).log.nLEDSepRemoved            = [sum(rem1), sum(rem2)]; % n samples removed per LED by separation check
obj.trialMetaData(trialIterator).log.PosLoadingStats(1,:)      = sum(~isnan(squeeze(ledPos(:,1,:))),1) / nSamp;

%%
% interpolate each LED across short gaps only - over short gaps both LEDs (and therefore position and head direction) are
% reliable. Longer gaps are left NaN here and dealt with on the level of the combined position (single LED fallback / position
% interpolation up to 'maxPosInterpolate'), where head direction stays NaN.
% NOTE: max gap is 1s. Tested against ground truth (real LED loss patterns imposed on a well tracked session): HD interpolated
% across 0.5-1s gaps is off by ~5deg (median), but ~4% of these samples are > 30deg off (95th prctile ~25deg); for gaps <= 0.5s
% the 95th prctile is <= 6deg. If you need very accurate HD sampling (e.g. HD cell tuning width) you might want to change this
% to 0.5s (costs HD coverage in sessions with poor tracking of one LED, e.g. ~6% of samples in r1010 novel_morph)
maxLEDGap = 1; % in s
for i = 1:2
    ledPos(:,:,i) = interpShortGaps(ledPos(:,:,i), round(maxLEDGap / frameInt));
end
%
obj.trialMetaData(trialIterator).log.PosLoadingStats(2,:)    = sum(~isnan(squeeze(ledPos(:,1,:))),1) / nSamp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xy, interpInd] = interpShortGaps(xy, maxGapSamp)
% linearly interpolate runs of NaNs (rows of xy) of up to maxGapSamp samples; runs at the start/end of the data are filled
% with the first/last valid sample. Longer runs stay NaN
interpInd = false(size(xy,1),1);
okInd     = find(~isnan(xy(:,1)));
if isempty(okInd) || numel(okInd) == size(xy,1); return; end
d         = diff([false; isnan(xy(:,1)); false]);
runStart  = find(d == 1);
runEnd    = find(d == -1) - 1;
for r = find((runEnd - runStart + 1) <= maxGapSamp)'
    interpInd(runStart(r):runEnd(r)) = true;
end
fillInd   = find(interpInd);
for j = 1:size(xy,2)
    xy(fillInd,j) = interp1(okInd, xy(okInd,j), fillInd, 'linear');
    xy(fillInd(fillInd > okInd(end)),j) = xy(okInd(end),j);
    xy(fillInd(fillInd < okInd(1)),j)   = xy(okInd(1),j);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [frameCount, sampleT] = fixFrameCounts(obj,trialIterator,frameCount,sampleT)
%%
% very rarely the frame counter (as well as the camera sample times) are corrupt from some time point onwards in a trial. That means from there onwards we cannot know anymore where potential missing frames occured - 
% if the n is low and your analysis doesn't require very high temporal accuracy just linearly interpolating these is prob. fine 
if any(diff(double(frameCount)) < 0)

    lastgoodInd = find(diff(double(frameCount)) > 1000,1,'first'); % it seems the corrupt samples have usually outlandish numbers (several orders of magnitude larger than normal frame count) 
    nMissFrames = sum(~ismember(1:frameCount(lastgoodInd),frameCount(1:lastgoodInd)));

    frameCount = [frameCount(1:lastgoodInd);(frameCount(lastgoodInd)+1:length(frameCount)+nMissFrames)'];

    if length(frameCount) < length(obj.spikeData.sampleT{trialIterator})
        % nMissFrames = sum(~ismember(1:frameCount(lastgoodInd),frameCount(1:lastgoodInd)));
        nFrameMissmatch = length(obj.spikeData.sampleT{trialIterator}) - (length(frameCount) + nMissFrames); 
        extraFrameInd   = round(linspace(double(frameCount(lastgoodInd+1)),length(frameCount)-1,nFrameMissmatch));

        for i = 1:length(extraFrameInd)
            frameCount(extraFrameInd(i):end) = frameCount(extraFrameInd(i):end) + 1; 
            sampleT(extraFrameInd(i):end)    = sampleT(extraFrameInd(i):end) + 1/obj.trialMetaData(trialIterator).posFs; 
        end
    else
        nFrameMissmatch = 0;
    end
    warning('scaNpix::npixUtils::loadPosNPix:FrameCount is corrupt from sample %i onwards. There are %i frames missing in remaining pos data. These were linearly interpolated - you should be aware of this!',frameCount(lastgoodInd),nFrameMissmatch);
    %
    obj.trialMetaData(trialIterator).log.frameCountCorruptFromSample = frameCount(lastgoodInd);
    obj.trialMetaData(trialIterator).log.nInterpSamplesCorruptFrames = nFrameMissmatch;
end
%%
% deal with missing syncs - we just treat them as missing frames - this is a bit of a headache as sometimes there are incomplete sync pulses at the point they drop off (so they miss in npix stream but not in pos stream). We need to deal
% with those
if ~isempty(obj.trialMetaData(trialIterator).missedSyncPulses)
    % first figure out if we have some extra frames in the pos stream
    missedSyncs                = obj.trialMetaData(trialIterator).missedSyncPulses; % [index of last good sync, n missing pulses, time of last good sync]
    posFs                      = obj.trialMetaData(trialIterator).posFs;
    [addPosFrames,posFrameInd] = deal(nan(1,size(missedSyncs,1)));
    totalNAddPosFrames         = 0;
    dSampleT                   = diff(sampleT);
    dFrameCount                = diff(double(frameCount));

    for i = 1:size(missedSyncs,1)
        % The camera is triggered by the sync pulses, so if they stopped the camera pauses as well, but without the frame count advancing (unlike dropped
        % frames). Sometimes the camera still catches a pulse or two that is missing in the npix stream, so we find the last frame before the camera
        % paused. We search a generous window around the gap, so we don't depend on exact alignment of camera and npix clocks (which drift apart).
        gapDur   = (missedSyncs(i,2) + 1) / posFs;
        candInd  = find( sampleT(1:end-1) > missedSyncs(i,3) - 1 & sampleT(1:end-1) < missedSyncs(i,3) + gapDur + 1 & dSampleT > 1.5/posFs & dFrameCount == 1 );
        if ~isempty(candInd)
            [~, maxInd]     = max(dSampleT(candInd));
            posFrameInd(i)  = candInd(maxInd); % last frame before camera paused
            addPosFrames(i) = double(frameCount(posFrameInd(i))) - missedSyncs(i,1) - totalNAddPosFrames; % frameCount is uint, which would clip negative values
        else
            % camera kept running, i.e. only the npix stream missed pulses - all frames are present in pos stream
            [~, posFrameInd(i)] = min(abs(sampleT - missedSyncs(i,3)));
            addPosFrames(i)     = missedSyncs(i,2);
        end
        if addPosFrames(i) < 0 || addPosFrames(i) > missedSyncs(i,2)
            warning('scaNpix::npixUtils::loadPosNPix:Couldn''t match pos frames to missing sync pulses (chunk %i: %i extra pos frames for %i missing syncs). Alignment of pos and spike data after %.1fs might be off - check this!', i, addPosFrames(i), missedSyncs(i,2), missedSyncs(i,3));
            addPosFrames(i) = min(max(addPosFrames(i), 0), missedSyncs(i,2));
        end
        totalNAddPosFrames = totalNAddPosFrames + addPosFrames(i);
    end
    % then update framecount accordingly
    for i = 1:size(obj.trialMetaData(trialIterator).missedSyncPulses,1)
        frameCount(posFrameInd(i)+1:end) = frameCount(posFrameInd(i)+1:end) + obj.trialMetaData(trialIterator).missedSyncPulses(i,2) - addPosFrames(i);  
    end
end

end

%         tmp = ledPos(~isnan(ledPos(:,1,i)),:,i);
% %     pathDists        = sqrt( diff(led(:,1,i),[],1).^2 + diff(led(:,2,i),[],1).^2 ) ./ ppm(1); % % distances in m
% %     tempSpeed        = pathDists ./ diff(sampleT); % m/s
% %     tempSpeed(end+1) = tempSpeed(end);
% %     speedInd = tempSpeed > obj.params('posMaxSpeed');
%         pathDists        = sqrt( diff(tmp(:,1),[],1).^2 + diff(tmp(:,2),[],1).^2 ) ./ ppm(1); % % distances in m
%         tempSpeed        = pathDists ./ mean(diff(sampleT)); %diff(sampleT(~isnan(ledPos(:,1,i)))); % m/s
%         tempSpeed(end+1) = tempSpeed(end);
%         speedInd = tempSpeed > maxSpeed;
%         tmp(speedInd,:) = NaN;
%         ledPos(~isnan(ledPos(:,1,i)),:,i) = tmp;

    % cs_NMissed = [0;cumsum(obj.trialMetaData(trialIterator).missedSyncPulses(:,2)-addPosFrames')];
    % missedSyncPosInd = obj.trialMetaData(trialIterator).missedSyncPulses(:,1) + cumsum(addPosFrames)' + cs_NMissed(1:end-1);
    % for i = 1:size(obj.trialMetaData(trialIterator).missedSyncPulses,1)
    %     frameCount(frameCount>missedSyncPosInd(i)) = frameCount(frameCount>missedSyncPosInd(i)) + obj.trialMetaData(trialIterator).missedSyncPulses(i,2) - addPosFrames(i);  
    % end

            % if sampleT(frameCount==obj.trialMetaData(trialIterator).missedSyncPulses(i,1)+1) - sampleT(frameCount==obj.trialMetaData(trialIterator).missedSyncPulses(i,1)) >= 1.1*1/obj.trialMetaData(trialIterator).posFs
        %     addPosFrames(i) = 0;
        % else
        %     addPosFrames(i) = find(diff(sampleT(find(frameCount<=obj.trialMetaData(trialIterator).missedSyncPulses(i,1)+totalNAddPosFrames):end))>1.5*1/obj.trialMetaData(trialIterator).posFs,1,'first') - 1;
        % end


         % missing_pos  = [1;find(diff(find(~isnan(currLED(:,1))))>1);length(currLED)];
    % missPosChunks = [missing_pos [missing_pos(2:end)+1;length(currLED)]];
    % pix = zeros(0,2);
    % for j = 1:size(missPosChunks,1)
    % 
    %     pix = [pix;currLED(missPosChunks(j,1):missPosChunks(j,2),:)];
    % end
    % 
    % 
    % missing_pos1  = find(isnan(currLED(:,1)));
    % chunkInd1      = diff(find([true,diff(missing_pos1')>1,true]));
    % % find those missing chunks where light was lost for too long (i.e. rat moved too far in between)
    % chunkInd      = diff(find([true,diff(missing_pos')>1,true]));
    % C             = mat2cell(missing_pos',1,chunkInd);
    % C(chunkInd1==1) = [];
    % C(cellfun(@(x) length(x)>10,C)|chunkInd1==1) = []; % remove single missing samples
    % 
    % pixMap = accumarray(currLED([C{:}],:),1,[max(currLED(:,1),[],'omitnan') max(currLED(:,2),[],'omitnan')]);
    % [xPix,yPix] = find(pixMap > 0.05*length([C{:}]));
    % badPixInd = ismember(currLED,[xPix,yPix],'rows');
    % currLED(badPixInd,:) = NaN;

    % remove positions that are flanked by NaNs - these are mostly dodgy and are spurious values that don't correspond to tracking the LEDs (we have to accept that we'll remove a few legit positions as well)
    % for j = 2
    %     remPosInd = 1;
    %     while ~isempty(remPosInd)
    %         trackedPosInd        = ~isnan(currLED(:,1));
    %         remPosInd            = find(conv(trackedPosInd,ones(2*j+1,1),'same') <= j & trackedPosInd);
    %         remPosInd(remPosInd < j+1 | remPosInd > length(currLED) - j) = [];
    %         %
    %         currLED(remPosInd,:) = NaN;
    %     end
    % end
