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

fprintf('Loading pos data for %s .......... ', obj.trialNames{trialIterator});


%% process data

% open pos data the format is [frame count, greenXY, redXY winSzX, winSzY, timeStamp possibly other Data ]
fName = dir(fullfile(obj.dataPath{trialIterator},'trackingData', '*.csv'));
if isempty(fName)
    warning(['scaNpix::loadPosNPix:Can''t find csv file in ' obj.dataPath{trialIterator} '. Come on mate.']);
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

% led data - same format as for dacq
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
    % sampleT          = (0:1/obj.trialMetaData(trialIterator).posFs:length(led)/obj.trialMetaData(trialIterator).posFs)';
    % sampleT          = sampleT(1:length(led)); % pretend we have perfect sampling
    sampleT    = sampleT + 1; % add 1 so interpolation won't fail - this data is essentally meaningless
    frameCount = [1;(length(led):-1:2)'+10e2]; % make a mock frame count that is corrupt from sample 1 onwards so we can use the fix in 'fixFrameCount' (in-line func.)
    obj.trialMetaData(trialIterator).BonsaiCorruptFlag = true;
    warning('scaNpix::loadPosNPix:Point Grey data corrupt!');
else
    frameCount       = csvData{1} - csvData{1}(1) + 1;
    obj.trialMetaData(trialIterator).BonsaiCorruptFlag = false;
end

% deal with missing syncs - we just treat them as missing frames - this is a bit of a headache as sometimes there are incomplete sync pulses at the point they drop off (so they miss in npix stream but not in pos stream). We need to deal
% with those
% if ~isempty(obj.trialMetaData(trialIterator).missedSyncPulses)
%     % first figure out if we have some extra frames in the pos stream
%     addPosFrames = nan(1,size(obj.trialMetaData(trialIterator).missedSyncPulses,1));
%     totalNAddPosFrames = 0;
%     for i = 1:size(obj.trialMetaData(trialIterator).missedSyncPulses,1)
%         if sampleT(frameCount==obj.trialMetaData(trialIterator).missedSyncPulses(i,1)+1) - sampleT(frameCount==obj.trialMetaData(trialIterator).missedSyncPulses(i,1)) >= 1.1*1/obj.trialMetaData(trialIterator).posFs
%             addPosFrames(i) = 0;
%         else
%             addPosFrames(i) = find(diff(sampleT(find(frameCount==obj.trialMetaData(trialIterator).missedSyncPulses(i,1)+totalNAddPosFrames):end))>1.5*1/obj.trialMetaData(trialIterator).posFs,1,'first') - 1;
%         end
%         totalNAddPosFrames = totalNAddPosFrames + addPosFrames(i);
%     end
%     % then update framecount accordingly
%     cs_NMissed = [0;cumsum(obj.trialMetaData(trialIterator).missedSyncPulses(:,2)-addPosFrames')];
%     missedSyncPosInd = obj.trialMetaData(trialIterator).missedSyncPulses(:,1) + cumsum(addPosFrames)' + cs_NMissed(1:end-1);
%     for i = 1:size(obj.trialMetaData(trialIterator).missedSyncPulses,1)
%         frameCount(frameCount>missedSyncPosInd(i)) = frameCount(frameCount>missedSyncPosInd(i)) + obj.trialMetaData(trialIterator).missedSyncPulses(i,2) - addPosFrames(i);  
%     end
% end

frameCount = fixFrameCounts(obj,trialIterator,frameCount,sampleT);

% deal with missing frames (if any) - this currently doesn't take into account if 1st frame(s) would be missing, but I am not sure this would
% actually ever happen (as 1st frame should always be triggered fine)
% first check if there are any...
missFrames       = find(~ismember(1:frameCount(end),frameCount));
nMissFrames      = length(missFrames);
if ~isempty(missFrames)
    fprintf('Note: There are %i missing frames in tracking data for %s.\n', nMissFrames, obj.trialMetaData(trialIterator).filename);
    
    temp                   = zeros(length(led)+nMissFrames, 2, obj.trialMetaData(trialIterator).nLEDs);
    temp(missFrames,:,:)   = nan;
    temp(temp(:,1)==0,:,:) = led;
    led                    = temp;
    
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
    ppm(:) = mean(envSzPix ./ (obj.trialMetaData(trialIterator).envSize ./ 100) );
else
    % this case should be default
    if ~circleFlag
        % recover all corner coords from 2 points - this should be independent of box misalignment with cam window
        knownDist = sqrt( (obj.trialMetaData(trialIterator).envBorderCoords(1,1)-obj.trialMetaData(trialIterator).envBorderCoords(1,2))^2 + (obj.trialMetaData(trialIterator).envBorderCoords(2,1)-obj.trialMetaData(trialIterator).envBorderCoords(2,2))^2 );
        ppm(:) = round( mean( knownDist ./ (sqrt(sum(obj.trialMetaData(trialIterator).envSize.^2)) ./ 100) ) );
        % full set
        obj.trialMetaData(trialIterator).envBorderCoords = scanpix.helpers.findBoxCorners(obj.trialMetaData(trialIterator).envBorderCoords(:,1),ppm(1)*(obj.trialMetaData(trialIterator).envSize(1)/100), obj.trialMetaData(trialIterator).envBorderCoords(:,2),ppm(1)*(obj.trialMetaData(trialIterator).envSize(2)/100));
        % now align env coords with the camera window
        for i = 1:2
            led(:,:,i) = scanpix.helpers.rotatePoints(led(:,:,i),[obj.trialMetaData(trialIterator).envBorderCoords(1,1),obj.trialMetaData(trialIterator).envBorderCoords(1,2);obj.trialMetaData(trialIterator).envBorderCoords(2,1),obj.trialMetaData(trialIterator).envBorderCoords(2,2)]);
        end
        %                     envSzPix  = [abs(obj.trialMetaData(trialIterator).envBorderCoords(1,1)-obj.trialMetaData(trialIterator).envBorderCoords(1,2)), abs(obj.trialMetaData(trialIterator).envBorderCoords(1,3)-obj.trialMetaData(trialIterator).envBorderCoords(2,3))];
    else
        [xCenter, yCenter, radius, ~] = scanpix.fxchange.circlefit(obj.trialMetaData(trialIterator).envBorderCoords(1,:), obj.trialMetaData(trialIterator).envBorderCoords(2,:));
        envSzPix = [2*radius 2*radius];
        ppm(:) = round( mean( envSzPix ./ (obj.trialMetaData(trialIterator).envSize ./ 100) ) );
    end
    %                 ppm(:) = round( mean( envSzPix ./ (obj.trialMetaData(trialIterator).envSize ./ 100) ) );
end

%% post process - basically as scanpix.dacqUtils.postprocess_data_v2
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
led = fixPositions(led, mean(diff(sampleT)), ppm(1), obj, trialIterator );

% smooth
kernel         = ones( ceil(obj.params('posSmooth') * obj.params('posFs')), 1)./ ceil( obj.params('posSmooth') * obj.params('posFs') ); % as per Ephys standard - 400ms boxcar filter
% Smooth lights individually, then get direction.
smLightFront   = imfilter(led(:, :, 1), kernel, 'replicate');
smLightBack    = imfilter(led(:, :, 2), kernel, 'replicate');

correction                              = obj.trialMetaData(trialIterator).LEDorientation(1); %To correct for light pos relative to rat subtract angle of large light
dirData                                 = mod((180/pi) * ( atan2(smLightFront(:,2)-smLightBack(:,2), smLightFront(:,1)-smLightBack(:,1)) ) - correction, 360); %
obj.posData(1).direction{trialIterator} = dirData(:);
% Get position from smoothed individual lights %%
wghtLightFront = 1-obj.params('posHead');
wghtLightBack  = obj.params('posHead');
xy = smLightFront .* wghtLightFront + smLightBack .* wghtLightBack;  %

% pos data
obj.posData(1).XYraw{trialIterator}        = xy;
obj.posData(1).XY{trialIterator}           = [floor(xy(:,1)) + 1, floor(xy(:,2)) + 1];
obj.posData(1).sampleT{trialIterator}      = sampleT; % this is redundant as we don't want to use the sample times from the PG camera

obj.trialMetaData(trialIterator).ppm       = ppm(1);
obj.trialMetaData(trialIterator).ppm_org   = ppm(2);

% scale position
boxExt = obj.trialMetaData(trialIterator).envSize / 100 * obj.trialMetaData(trialIterator).ppm;
scanpix.maps.scalePosition(obj, trialIterator,'envszpix', boxExt,'circflag',circleFlag); % need to enable this for circular env as well!

% running speed
pathDists                                  = sqrt( diff(xy(:,1)).^2 + diff(xy(:,2)).^2 ) ./ ppm(1) .* 100; % distances in cm
obj.posData(1).speed{trialIterator}        = pathDists ./ diff(sampleT); % cm/s
obj.posData(1).speed{trialIterator}(end+1) = obj.posData(1).speed{trialIterator}(end);

fprintf('  DONE!\n');

end


function ledPos = fixPositions(ledPos,sampleT,ppm,obj,trialIterator)

obj.trialMetaData(trialIterator).log.PosLoadingStats(1,:) = sum(~isnan(squeeze(ledPos(:,1,:))),1) / size(ledPos,1);

% first we remove positions that are flanked by NaNs - these are mostly dodgy and are spurious values that don't correspond to tracking the LEDs (we have to accept that we'll remove a few legit positions as well)
for i = 1:2
    currLED = ledPos(:,:,i);
    remPosInd = 1;
    while ~isempty(remPosInd) %any(speedInd)
        trackedPosInd        = ~isnan(currLED(:,1));
        remPosInd            = find(conv(trackedPosInd,ones(5,1),'same') <= 2 & trackedPosInd);
        currLED(remPosInd,:) = NaN;
    end  
    % now look for tracking errors by speed - we'll ignore all the NaNs here as these prevent to identify some dodgy samples (again we might lose a few legit samples here when the light wasn't tracked for too
    % long continuously)
    validPos                        = currLED(~isnan(currLED(:,1)),:);
    pathDists                       = sqrt( diff(validPos(:,1),[],1).^2 + diff(validPos(:,2),[],1).^2 ) ./ ppm(1); % % distances in m
    tempSpeed                       = pathDists ./ sampleT; %diff(sampleT(~isnan(ledPos(:,1,i)))); % m/s
    tempSpeed(end+1)                = tempSpeed(end);
    speedInd                        = tempSpeed > obj.params('posMaxSpeed');
    validPos(speedInd,:)            = NaN;
    currLED(~isnan(currLED(:,1)),:) = validPos;
    ledPos(:,:,i)                   = currLED;
end

% interpolate between good samples  
for i = 1:2
    % find all missing positions/led
    missing_pos   = find(isnan(ledPos(:,1,i)));
    % find those missing chunks where light was lost for too long (i.e. rat moved too far in between)
    idx           = find(diff(missing_pos)>1);
    idx           = idx(1:end-1); % can ignore last entry
    missPosChunks = [[max([1,missing_pos(1)-1]); missing_pos(idx(1:end-1)+1)-1],missing_pos(idx)+1]; % make sure first index~=0
    indTooLong    = sqrt(diff([ledPos(missPosChunks(:,1),1,i),ledPos(missPosChunks(:,2),1,i)],[],2).^2+diff([ledPos(missPosChunks(:,1),2,i),ledPos(missPosChunks(:,2),2,i)],[],2).^2) ./ ppm .* 100 > obj.params('maxPosInterpolate');
    missPosChunks = missPosChunks(indTooLong,:); % only keep these
    % remove all bad chunks
    for j = 1:size(missPosChunks,1)
        missing_pos = missing_pos(~ismember(missing_pos,missPosChunks(j,1)+1:missPosChunks(j,2)-1));
    end
    % interpolate as per usual
    ok_pos      = find(~isnan(ledPos(:,1,i)));
    for j = 1:2
        ledPos(missing_pos, j, i)                            = interp1(ok_pos, ledPos(ok_pos, j, i), missing_pos, 'linear');
        ledPos(missing_pos(missing_pos > max(ok_pos)), j, i) = ledPos( max(ok_pos), j, i);
        ledPos(missing_pos(missing_pos < min(ok_pos)), j, i) = ledPos( min(ok_pos), j, i);
    end
    obj.trialMetaData(trialIterator).log.PosLoadingStats(2,i) = (length(missing_pos)+length(ok_pos)) / size(ledPos,1);
end
% now last sanity check - check for samples where distance between LEDs is too large - set these position for the LED that was tracked worse overall to the ones from the better tracked light 
[~, maxInd] = max([sum(~isnan(ledPos(:,1,1))),sum(~isnan(ledPos(:,1,2)))]);
LEDdistInd = sqrt( (ledPos(:,1,1) - ledPos(:,1,2)).^2 + (ledPos(:,2,1) - ledPos(:,2,2)).^2 ) ./ ppm(1) .* 100 > 15; % 
ledPos(LEDdistInd,:,maxInd~=[1 2]) = ledPos(LEDdistInd,:,maxInd); 

end

function frameCount = fixFrameCounts(obj,trialIterator,frameCount,sampleT)

% very rarely the frame counter (as well as the camera sample times) are corrupt from some time point onwards in a trial. That means from there onwards we cannot know anymore where potential missing frames occured - 
% if the n is low and your analysis doesn't require very high temporal accuracy just linearly interpolating these is prob. fine 
if any(diff(double(frameCount)) < 0)

    lastgoodInd = find(diff(double(frameCount)) > 1000,1,'first'); % it seems the corrupt samples have ussually outlandish numbers (several orders of magnitude larger than normal frame count) 
    nMissFrames = sum(~ismember(1:frameCount(lastgoodInd),frameCount(1:lastgoodInd)));

    frameCount = [frameCount(1:lastgoodInd);(frameCount(lastgoodInd)+1:length(frameCount)+nMissFrames)'];

    if length(frameCount) < length(obj.spikeData.sampleT{trialIterator})
        % nMissFrames = sum(~ismember(1:frameCount(lastgoodInd),frameCount(1:lastgoodInd)));
        nFrameMissmatch = length(obj.spikeData.sampleT{trialIterator}) - (length(frameCount) + nMissFrames); 
        extraFrameInd = round(linspace(double(frameCount(lastgoodInd+1)),length(frameCount)-1,nFrameMissmatch));

        for i = 1:length(extraFrameInd)
            frameCount(extraFrameInd(i):end) = frameCount(extraFrameInd(i):end) + 1; 
        end
    else
        nFrameMissmatch = 0;
    end
    warning('scaNpix::loadPosNPix:FrameCount is corrupt from sample %i onwards. There are %i frames missing in remaining pos data. These were linearly interpolated - you should be aware of this!',frameCount(lastgoodInd),nFrameMissmatch);
    %
    obj.trialMetaData(trialIterator).log.frameCountCorruptFromSample = frameCount(lastgoodInd);
    obj.trialMetaData(trialIterator).log.nInterpSamplesCorruptFrames = nFrameMissmatch;
end


% deal with missing syncs - we just treat them as missing frames - this is a bit of a headache as sometimes there are incomplete sync pulses at the point they drop off (so they miss in npix stream but not in pos stream). We need to deal
% with those
if ~isempty(obj.trialMetaData(trialIterator).missedSyncPulses)
    % first figure out if we have some extra frames in the pos stream
    [addPosFrames,posFrameInd] = deal(nan(1,size(obj.trialMetaData(trialIterator).missedSyncPulses,1)));
    totalNAddPosFrames = 0;

    for i = 1:size(obj.trialMetaData(trialIterator).missedSyncPulses,1)
        % find the relevant pos frame where the syncs are missing - 'min' should be fine here as next sampleT will correspond to time when syncs came back, so there should be a temporal gap
        posFrameInd(i) = find(abs(sampleT - obj.trialMetaData(trialIterator).missedSyncPulses(i,3)) < 1/obj.trialMetaData(trialIterator).posFs,1,'last');
        %
        addPosFrames(i) = frameCount(posFrameInd(i)) - obj.trialMetaData(trialIterator).missedSyncPulses(i,1) - totalNAddPosFrames;
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
