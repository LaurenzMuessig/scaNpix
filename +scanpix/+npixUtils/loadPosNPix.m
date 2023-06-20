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
    disp(['Can''t find csv file in ' obj.dataPath{trialIterator} '. Come on mate.']);
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
    sampleT          = (0:1/obj.trialMetaData(trialIterator).posFs:length(led)/obj.trialMetaData(trialIterator).posFs)';
    sampleT          = sampleT(1:length(led)); % pretend we have perfect sampling
    frameCount       = 1:length(led); % pretend we are not missing any frames
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
    
    
%     nMissedPulses = floor((syncTTLs(missedSyncs+1) - syncTTLs(missedSyncs)) * obj.params('posFs'));
%     missedPulses  = missedSyncs+1:missedSyncs+nMissedPulses;
%     interp_pulseT          = interp1([1:missedSyncs,missedSyncs+nMissedPulses+1:length(syncTTLs)+nMissedPulses], syncTTLs', missedPulses);
%     temp                   = zeros(length(syncTTLs)+nMissedPulses,1);
%     temp(missedPulses,1)   = interp_pulseT;
%     temp(temp(:,1) == 0,1) = syncTTLs;
%     syncTTLs               = temp;
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
    obj.trialMetaData(trialIterator).objectPos = obj.trialMetaData(trialIterator).objectPos .* scaleFact;
    obj.trialMetaData(trialIterator).envBorderCoords = obj.trialMetaData(trialIterator).envBorderCoords .* scaleFact;
    if circleFlag
        [xCenter, yCenter, radius] = deal(xCenter*scaleFact,yCenter*scaleFact,radius*scaleFact);
    end
end

% remove tracking errors (i.e. too fast)
for i = 1:2
    % speed
    pathDists        = sqrt( diff(led(:,1,i),[],1).^2 + diff(led(:,2,i),[],1).^2 ) ./ ppm(1); % % distances in m
    tempSpeed        = pathDists ./ diff(sampleT); % m/s
    tempSpeed(end+1) = tempSpeed(end);
    speedInd = tempSpeed > obj.params('posMaxSpeed');
    % env borders
    if ~circleFlag
        envSzInd = led(:,1,i) < 0.95 * min(obj.trialMetaData(trialIterator).envBorderCoords(1,:)) | led(:,1,i) > 1.05 * max(obj.trialMetaData(trialIterator).envBorderCoords(1,:)) | led(:,2,i) < 0.95 * min(obj.trialMetaData(trialIterator).envBorderCoords(2,:)) | led(:,2,i) > 1.05 * max(obj.trialMetaData(trialIterator).envBorderCoords(2,:));
    else
        envSzInd = (led(:,1,i) - xCenter).^2 + (led(:,2,i) - yCenter).^2 > radius^2; % points outside of environment
    end
    % filter out
    led(speedInd | envSzInd,:,i) = NaN;
end

% interpolate missing positions
for i = 1:2
    missing_pos = find(isnan(led(:,1,i)));
    ok_pos      = find(~isnan(led(:,1,i)));
    for j = 1:2
        led(missing_pos, j, i)                            = interp1(ok_pos, led(ok_pos, j, i), missing_pos, 'linear');
        led(missing_pos(missing_pos > max(ok_pos)), j, i) = led( max(ok_pos), j, i);
        led(missing_pos(missing_pos < min(ok_pos)), j, i) = led( min(ok_pos), j, i);
    end
end

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
xy = (smLightFront .* wghtLightFront + smLightBack .* wghtLightBack);  %

% pos data
obj.posData(1).XYraw{trialIterator}        = xy;
obj.posData(1).XY{trialIterator}           = [double( floor(xy(:,1)) + 1 ), double( floor(xy(:,2)) + 1 )];
obj.posData(1).sampleT{trialIterator}      = sampleT; % this is redundant as we don't want to use the sample times from the PG camera

obj.trialMetaData(trialIterator).ppm       = ppm(1);
obj.trialMetaData(trialIterator).ppm_org   = ppm(2);

% scale position
boxExt = obj.trialMetaData(trialIterator).envSize / 100 * obj.trialMetaData(trialIterator).ppm;
scanpix.maps.scalePosition(obj, trialIterator,'envszpix', boxExt,'circflag',circleFlag); % need to enable this for circular env as well!

% running speed
%             pathDists                                  = sqrt( (obj.posData(1).XY{trialIterator}(1:end-1,1) - obj.posData(1).XY{trialIterator}(2:end,1)).^2 + (obj.posData(1).XY{trialIterator}(1:end-1,2) - obj.posData(1).XY{trialIterator}(2:end,2)).^2 ) ./ ppm(1) .* 100; % distances in cm
pathDists                                  = sqrt( diff(xy(:,1)).^2 + diff(xy(:,2)).^2 ) ./ ppm(1) .* 100; % distances in cm
obj.posData(1).speed{trialIterator}        = pathDists ./ diff(sampleT); % cm/s
obj.posData(1).speed{trialIterator}(end+1) = obj.posData(1).speed{trialIterator}(end);

fprintf('  DONE!\n');

end

