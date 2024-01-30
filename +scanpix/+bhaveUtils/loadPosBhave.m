function loadPosBhave(obj, trialIterator)
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

% open pos data the format is [frame count,timeStamp, posX, posY, area, possibly other Data ]
fID = fopen(fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} obj.fileType]),'rt');
header = textscan(fID,'%s',1);
nColumns = length(strsplit(header{1}{1},','));
fmt = '%u%f%f%f%f%f%u%u';
% allow for any n of additonal fields from Bonsai output
if nColumns > 8; fmt = [fmt repmat('%u',1,nColumns-8)]; end

csvData = textscan(fID,fmt,'HeaderLines',1,'delimiter',',');
fclose(fID);

% add possible extra data from Bonsai
if nColumns > 8
    bhaveData = csvData(9:end);
end

% 
pos         = [csvData{3}, csvData{4}];
pos(pos==0) = NaN;

% sample Times
timeStamps  = csvData{2};
sampleT     = scanpix.npixUtils.convertPointGreyCamTimeStamps(timeStamps); % starts @ 0

% in case logging point grey data was corrupt
if all(sampleT == 0)
    sampleT          = (0:1/obj.trialMetaData(trialIterator).posFs:length(pos)/obj.trialMetaData(trialIterator).posFs)';
    sampleT          = sampleT(1:length(pos)); % pretend we have perfect sampling
    frameCount       = 1:length(pos); % pretend we are not missing any frames
    obj.trialMetaData(trialIterator).BonsaiCorruptFlag = true;
    warning('Point Grey data corrupt!');
else
    frameCount       = csvData{1} - csvData{1}(1) + 1;
    obj.trialMetaData(trialIterator).BonsaiCorruptFlag = false;
end

% deal with missing frames (if any) - this currently doesn't take into
% account if 1st frame(s) would be missing, but I am not sure this would
% actually ever happen (as 1st frame should always be triggered fine)
% first check if there are any...
missFrames       = find(~ismember(1:frameCount(end),frameCount));
nMissFrames      = length(missFrames);
if ~isempty(missFrames)
    fprintf('Note: There are %i missing frames in tracking data for %s.\n', nMissFrames, obj.trialMetaData(trialIterator).filename);
    
    temp                   = zeros(length(pos)+nMissFrames, 2);
    temp(missFrames,:,:)   = nan;
    temp(temp(:,1)==0,:,:) = pos;
    pos                    = temp;
    
    % interpolate sample times
    interp_sampleT           = interp1(double(frameCount), sampleT, missFrames);
    temp2                    = zeros(length(pos),1);
    temp2(missFrames,1)      = interp_sampleT;
    temp2(temp2(:,1) == 0,1) = sampleT;
    sampleT                  = temp2;
    
    % add mising trials to behav data
    if nColumns > 8
        for i = 1:size(bhaveData,2)
            temp                 = zeros(length(pos), size(bhaveData{i},2));
            temp(missFrames,:)   = nan;
            temp(temp(:,1)==0,:) = bhaveData{i};
            bhaveData{i}         = temp;
        end
    end

 
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
        pos = scanpix.helpers.rotatePoints(pos,[obj.trialMetaData(trialIterator).envBorderCoords(1,1),obj.trialMetaData(trialIterator).envBorderCoords(1,2);obj.trialMetaData(trialIterator).envBorderCoords(2,1),obj.trialMetaData(trialIterator).envBorderCoords(2,2)]);
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
    pos = floor(pos .* scaleFact);
    ppm(1) = obj.params('ScalePos2PPM');
    obj.trialMetaData(trialIterator).envBorderCoords = obj.trialMetaData(trialIterator).envBorderCoords .* scaleFact;
    if circleFlag
        [xCenter, yCenter, radius] = deal(xCenter*scaleFact, yCenter*scaleFact, radius*scaleFact);
    end
    obj.trialMetaData(trialIterator).PosIsScaled = true;
else
    obj.trialMetaData(trialIterator).PosIsScaled = false;
end

% remove tracking errors that fall outside box
% env borders
borderTolerancePix = ppm(1)/100*2.5; % we'll assume 1 standard rate map bin tolerance
if ~circleFlag
    envSzInd = pos(:,1) < min(obj.trialMetaData(trialIterator).envBorderCoords(1,:))-borderTolerancePix | pos(:,1) > max(obj.trialMetaData(trialIterator).envBorderCoords(1,:))+borderTolerancePix | pos(:,2) < min(obj.trialMetaData(trialIterator).envBorderCoords(2,:))-borderTolerancePix | pos(:,2) > max(obj.trialMetaData(trialIterator).envBorderCoords(2,:))+borderTolerancePix;
else
    envSzInd = (pos(:,1) - xCenter).^2 + (pos(:,2) - yCenter).^2 > (radius+borderTolerancePix)^2; % points outside of environment
end
% filter out
pos(envSzInd,:) = NaN;

% fix positions (inline subfunction)
pos = fixPositions(pos, mean(diff(sampleT)), ppm(1), obj.params('posMaxSpeed'), obj.params('maxPosInterpolate'));

% smooth
kernel = ones( ceil(obj.params('posSmooth') * obj.params('posFs')), 1)./ ceil( obj.params('posSmooth') * obj.params('posFs') ); % as per Ephys standard - 400ms boxcar filter
% Smooth lights individually, then get direction.
xy = imfilter(pos, kernel, 'replicate');
% movement direction!
dirData                              = mod((180/pi) * ( atan2(xy(2:end,2)-xy(1:end-1,2), xy(2:end,1)-xy(1:end-1,1)) ), 360); % not sure this is right
dirData(end+1)                       = dirData(end); 
obj.posData.direction{trialIterator} = dirData(:);

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

% crop overhang at the end
endIdxNPix                           = min( [ length(obj.posData.sampleT{trialIterator}), find(obj.posData.sampleT{trialIterator} < obj.trialMetaData(trialIterator).duration,1,'last') + 1]);
obj.posData.XYraw{trialIterator}     = obj.posData.XYraw{trialIterator}(1:endIdxNPix,:);
obj.posData.XY{trialIterator}        = obj.posData.XY{trialIterator}(1:endIdxNPix,:);
obj.posData.speed{trialIterator}     = obj.posData.speed{trialIterator}(1:endIdxNPix,:);
obj.posData.direction{trialIterator} = obj.posData.direction{trialIterator}(1:endIdxNPix,:);
obj.posData.sampleT{trialIterator}   = obj.posData.sampleT{trialIterator}(1:endIdxNPix,:);

% add possible extra data from Bonsai
if nColumns > 8
    obj.bhaveData(1).data(trialIterator) = cellfun(@(x) x(1:endIdxNPix,:),csvData(9:end),'uni',0);
end

fprintf('  DONE!\n');

end


function pos = fixPositions(pos,sampleT,ppm,maxSpeed,maxPosInterpolateCM)

% first we remove positions that are flanked by NaNs - these are mostly dodgy and are spurious values that don't correspond to tracking the LEDs (we have to accept that we'll remove a few legit positions as well)
remPosInd = 1;
while ~isempty(remPosInd) %any(speedInd)
    trackedPosInd = ~isnan(pos(:,1));
    remPosInd = find(conv(trackedPosInd,ones(5,1),'same') <= 2 & trackedPosInd);
    pos(remPosInd,:) = NaN;
end
% now look for tracking errors by speed - we'll ignore all the NaNs here as these prevent to identify some dodgy samples (again we might lose a few legit samples here when the light wasn't tracked for too
% long continuously)
validPos                = pos(~isnan(pos(:,1)),:);
pathDists               = sqrt( diff(validPos(:,1),[],1).^2 + diff(validPos(:,2),[],1).^2 ) ./ ppm(1); % % distances in m
tempSpeed               = pathDists ./ sampleT; %diff(sampleT(~isnan(ledPos(:,1,i)))); % m/s
tempSpeed(end+1)        = tempSpeed(end);
speedInd                = tempSpeed > maxSpeed;
validPos(speedInd,:)    = NaN;
pos(~isnan(pos(:,1)),:) = validPos;


% interpolate between good samples  
% find all missing positions/led
missing_pos   = find(isnan(pos(:,1)));
if ~isempty(missing_pos) && length(missing_pos) > 1
    % find those missing chunks where light was lost for too long (i.e. rat moved too far in between)
    
%     idx           = find(diff(missing_pos)>1);
    idx = diff(missing_pos);
    if sum(idx) - max(idx) == length(idx) - 1
        missPosChunks = [missing_pos(find(idx==1,1,'first')) missing_pos(find(idx==1,1,'last'))]; % only 1 valid chunk
    else
        idx = find(idx(1:end-1)>1);
        missPosChunks = [[max([1,missing_pos(1)-1]); missing_pos(idx(1:end-1)+1)-1],missing_pos(idx)+1]; % make sure first index~=0
    end
%     missPosChunks = [[max([1,missing_pos(1)-1]); missing_pos(idx(1:end-1)+1)-1],missing_pos(idx)+1]; % make sure first index~=0
    if ~isempty(missPosChunks)
        indTooLong    = sqrt(diff([pos(missPosChunks(:,1),1),pos(missPosChunks(:,2),1)],[],2).^2+diff([pos(missPosChunks(:,1),2),pos(missPosChunks(:,2),2)],[],2).^2) ./ ppm .* 100 > maxPosInterpolateCM;
        missPosChunks = missPosChunks(indTooLong,:); % only keep these
        % remove all bad chunks
        for i = 1:size(missPosChunks,1)
            missing_pos = missing_pos(~ismember(missing_pos,missPosChunks(i,1)+1:missPosChunks(i,2)-1));
        end
    end
end

% interpolate as per usual
ok_pos      = find(~isnan(pos(:,1)));
for i = 1:2
    pos(missing_pos, i)                            = interp1(ok_pos, pos(ok_pos, i), missing_pos, 'linear');
    pos(missing_pos(missing_pos > max(ok_pos)), i) = pos( max(ok_pos), i);
    pos(missing_pos(missing_pos < min(ok_pos)), i) = pos( min(ok_pos), i);
end

end

