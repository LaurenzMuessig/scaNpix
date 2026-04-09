function loadPos(obj, trialIterator)
% loadPos - load position data from DACQ files
%
% Syntax:  loadPos(obj, trialIterator)
%
% Inputs:
%    obj           - ephys class object ('dacq')
%    trialIterator - numeric index for trial to be loaded
%
% Outputs:
%
% See also:
%
% TW/LM 2020 (adapted from org. SCAN function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialIterator (1,1) {mustBeNumeric}
end

%%
fprintf('Loading pos file for %s [..........] ', obj.trialNames{trialIterator});

%%% Read data from File %%%
fid = fopen(fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.pos']),'r','ieee-be');
if fid ==-1
    ME = MException('scaNpix::load_pos:posFileNotFound', ['Could not open pos file ''' fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.pos']) '']);
    throw(ME);
end

C = textscan(fid,'%s %[^\r\n]', 27);  % CAUTION! N Lines of header hard-coded.
frewind(fid);
posHeader    = [cat(1,C{1}) cat(1,C{2})];
% Look for data start marker
headerText   = fread(fid,800,'int8');
headerText   = char(abs(headerText))';
ds           = strfind(headerText,'data_start') + 10;
% Read data %
fseek(fid,(ds-1),'bof');                                    % Set file position indicator at data start
nPosSamples  = sscanf(scanpix.dacqUtils.getValue(posHeader,'num_pos_samples'), '%d');
data         = fread(fid, nPosSamples * 10, 'uint16');
fclose(fid);

%%% Reshape into correct output format %%%
data         = reshape(data, [10 nPosSamples])';
data         = data(:,3:10);
data         = reshape(data, [nPosSamples, 2, 4]); % Now in format [nSamp, x:y, nLight]
% Separate numpix, if existing, switch format to [nSamp, nLight, x:y] %
nLights      = obj.trialMetaData(trialIterator).tracked_spots;
if strcmp(scanpix.dacqUtils.getValue(posHeader,'pos_format'), 't,x1,y1,x2,y2,numpix1,numpix2') && nLights <= 2
    led_pos  = nan(nPosSamples, nLights, 2);
    led_pix  = nan(nPosSamples, nLights);
    for i = 1:nLights
        for j = 1:2
            led_pos(:,i,j) = data(:,j,i);
        end
        led_pix(:,i) = data(:,i,3); % numpix always seems to start at 3rd light (5th entry)
    end
else
    % not sure this is still necessary as we are not allowing legacy data
    % anyway
    led_pos = nan(nPosSamples, nLights, 2);
    led_pix = [];
    for i = 1:nLights
        for j = 1:2
            led_pos(:,i,j) = data(:,j,i);
        end
    end
end
led_pos(led_pos==1023) = NaN;   % mTint functions are
led_pix(led_pos==1023) = NaN;   % expecting this.

% scale pos data to specific PPM
obj.trialMetaData(trialIterator).ppm_org = sscanf(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'),'%d');

if ~isempty(obj.params('ScalePos2PPM'))
    scaleFact = (obj.params('ScalePos2PPM') / str2double(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre')));
    led_pos = floor( led_pos .*  scaleFact);
    % posHeader{strcmp(posHeader(:,1),'pixels_per_metre'),2} = num2str(obj.params('ScalePos2PPM'));
    obj.trialMetaData(trialIterator).ppm = obj.params('ScalePos2PPM');

    %
    obj.trialMetaData(trialIterator).PosIsScaled = true;
else
    obj.trialMetaData(trialIterator).PosIsScaled = false;
    obj.trialMetaData(trialIterator).ppm         = obj.trialMetaData(trialIterator).ppm_org;
end
% obj.trialMetaData(trialIterator).ppm = sscanf(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'),'%d');

% post process
obj.posData(1).XYraw{trialIterator}    = led_pos;
obj.params('posFs')                    = sscanf(scanpix.dacqUtils.getValue(posHeader,'sample_rate'),'%d');
obj.trialMetaData(trialIterator).posFs = sscanf(scanpix.dacqUtils.getValue(posHeader,'sample_rate'),'%d');
%Interpolate, Filter and Smooth %  %%%%%%%%%% THIS COULD USE A BIT MORE EFFICIENCY AND BE UPDATED
% [obj.posData(1).XYraw{trialIterator}, obj.posData(1).direction{trialIterator}, obj.posData(1).speed{trialIterator}, obj.trialMetaData(trialIterator).log.PosLoadingStats] = scanpix.dacqUtils.postprocess_pos_data_v2(posFile, obj.params('posMaxSpeed'), obj.params('posSmooth'), obj.trialMetaData(trialIterator), obj.params('posHead'));
postprocess_posData(obj,trialIterator,led_pos,led_pix);

% remove DACQ overhang
if obj.trialMetaData(trialIterator).duration * obj.params('posFs') < length(led_pos)
    obj.posData(1).XYraw{trialIterator}     = obj.posData(1).XYraw{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs'),:); % truncate data
    obj.posData(1).XY{trialIterator}        = obj.posData(1).XY{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs'),:); % truncate data
    obj.posData(1).direction{trialIterator} = obj.posData(1).direction{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs')); % truncate data
    obj.posData(1).speed{trialIterator}     = obj.posData(1).speed{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs')); % truncate data
end

% convert to integers
% obj.posData(1).XY{trialIterator} = [double( floor(obj.posData(1).XYraw{trialIterator}(:,1)) + 1 ), double( floor(obj.posData(1).XYraw{trialIterator}(:,2)) + 1 )];  %%% NECESSARY??

%
if ~obj.params('scalePos2CamWin') && ~isempty(obj.trialMetaData(trialIterator).envSize )
    boxExt = obj.trialMetaData(trialIterator).envSize / 100 * obj.trialMetaData(trialIterator).ppm;
    scanpix.maps.scalePosition(obj, trialIterator,'envszpix', boxExt);
    %
    % obj.posData(1).XY{trialIterator} = scanpix.helpers.rotatePoints(obj.posData(1).XY{trialIterator});
end

fprintf('  DONE!\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  postprocess_posData(obj,trialIterator,led_pos,led_pix )

%%
% stats = nan(3,2);
n_pos  = size(led_pos,1); % Don't use led_pix because this is absent from the older format
n_leds = size(led_pos,2);

% pos_sample_rate = sscanf(scanpix.dacqUtils.getValue(posdata.header, 'sample_rate'),'%d');

% For 2 spot tracking, check for instances of swapping (often happens if one LED bright, one less bright).
if n_leds == 2 && ~(isempty(led_pix)) % Only check if we actually have led_pix
    swap_list                        = led_swap_filter(led_pos, led_pix);
    dum                              = led_pos(swap_list, 1, :);
    led_pos(swap_list, 1, :) = led_pos(swap_list, 2, :);
    led_pos(swap_list, 2, :) = dum;
    % dum                              = led_pix(swap_list, 1);
    % led_pix(swap_list, 1)    = led_pix(swap_list, 2);
    % led_pix(swap_list, 2)    = dum;
end
%
% stats(3,1) = length(swap_list) / size(led_pos,1);

% Filter points for those moving impossibly fast and set led_pos_x to NaN
max_pix_per_sample = obj.params('posMaxSpeed')*obj.trialMetaData(trialIterator).ppm/obj.trialMetaData(trialIterator).posFs ;
for led = 1:n_leds
    [n_jumpy, led_pos] = led_speed_filter( led_pos, max_pix_per_sample, led);
    if( n_jumpy > n_pos/3 )
        warning('scaNpix::dacqUtils::postprocess_posData: %d/%d positions rejected for led %d\n', n_jumpy, n_pos, led);
    end
end
%
obj.trialMetaData(trialIterator).log.PosLoadingStats(1,:) = sum(~isnan(squeeze(led_pos(:,:,1))),1) / size(led_pos,1);

% Interpolate to replace all NaN led positions
% 1/12/09 AJ: I've made changes to the following lines to make the pos.xy
% output more robust. interp1 ignores NaNs at the endpoints of a vector. So
% two lines have been added to take care of these. missing_pos is formed to
% reject poses where the system has failed to record either x OR y
% coordinate and ok_pos requires both x AND y coords to be present.
for led = 1:n_leds
        missing_pos = find(isnan(led_pos(:, led, 1)) | isnan(led_pos(:, led, 2)));
        ok_pos      = find(~isnan(led_pos(:, led, 1)) & ~isnan(led_pos(:, led, 2)));
    for k = 1:2
        led_pos(missing_pos, led, k)                          = interp1(ok_pos, led_pos(ok_pos, led, k), missing_pos, 'linear');
        led_pos(missing_pos(missing_pos>max(ok_pos)), led, k) = led_pos(max(ok_pos), led, k);
        led_pos(missing_pos(missing_pos<min(ok_pos)), led, k) = led_pos(min(ok_pos), led, k);
    end
end
%
obj.trialMetaData(trialIterator).log.PosLoadingStats(2,:) = sum(~isnan(squeeze(led_pos(:,:,1))),1) / size(led_pos,1);

% Estimate position, direction and speed from led_pos, using smoothing
kernel = ones(ceil( obj.params('posSmooth')*obj.trialMetaData(trialIterator).posFs), 1)./ceil( obj.params('posSmooth')*obj.trialMetaData(trialIterator).posFs);

% Need to know angles (and distances) of LEDs from rat (in bearing_colour_i in .pos header)
pos  = 1:n_pos;
pos2 = 1:n_pos-1;
xy   = nan(n_pos, 2);
if n_leds == 1
    xy(pos, :) = led_pos(pos, 1, :);
    xy         = imfilter(xy, kernel, 'replicate');
    % y from tracker increases downwards. dir is positive anticlockwise from X axis, 0<= dir <360
    %CB verified calc of dir in next line is correct & concords with tint
    dir(pos2)  = mod((180/pi)*(atan2( -xy(pos2+1, 2)+xy(pos2, 2), xy(pos2+1, 1)-xy(pos2, 1) )), 360);
    dir(n_pos) = dir(n_pos-1)';
    obj.posData(1).direction{trialIterator} = dir(:);
    % dir_disp   = dir'; %Return dir_disp for completness even though == dir
elseif n_leds == 2
    % 2 LEDs are assumed to be at 180deg to each other with their midpoint over the
    % animals head. Convention is that colbearings_set(1) is the large light (normally
    % at front) and (2) is the small light. Position of lights (defined in set
    % file header) is defined in degs with 0 being straight ahead of animal, values
    % increase anti-clockwise.
    
    % Smooth lights individually, then get direction. %% TW, 12/09/08: This method replicates TINT.
    smLightFront(pos, :) = imfilter(led_pos(pos, 1, :), kernel, 'replicate');
    smLightBack(pos, :)  = imfilter(led_pos(pos, 2, :), kernel, 'replicate');
    
    dir        = mod((180/pi)*( atan2(-smLightFront(pos,2)+smLightBack(pos,2), +smLightFront(pos,1)-smLightBack(pos,1)) ) - obj.trialMetaData(trialIterator).lightBearing(1), 360);
    obj.posData(1).direction{trialIterator} = dir(:);
    % Get position from smoothed individual lights %%  % TW, 12/09/08
    wghtLightFront = 1- obj.params('posHead');
    wghtLightBack  = obj.params('posHead');
    xy(pos, :)     = (smLightFront(pos, :).*wghtLightFront + smLightBack(pos, :).*wghtLightBack);  %%% CB added code todo headPos other than 0.5.
    obj.posData(1).XY{trialIterator} = xy;

    %Get heading from displacement too
    % dir_disp(pos2)  = mod((180/pi)*(atan2( -xy(pos2+1, 2)+xy(pos2, 2), xy(pos2+1, 1)-xy(pos2, 1) )), 360);
    % dir_disp(n_pos) = dir_disp(n_pos-1);
    % dir_disp        = dir_disp(:);
end

%%% Calculate speed, based on distance(sampleN+1-sampleN) %%%
speed(pos2)  = sqrt((xy(pos2+1,1)-xy(pos2,1)).^2+(xy(pos2+1,2)-xy(pos2,2)).^2);
speed(n_pos) = speed(n_pos-1);
speed        = speed.*(100*obj.trialMetaData(trialIterator).posFs/obj.trialMetaData(trialIterator).ppm);
obj.posData(1).speed{trialIterator} = speed(:);

% times = (1:n_pos)/pos_sample_rate;
% times = times(:);
end


function [n_jumpy, led_pos] = led_speed_filter(led_pos, max_pix_per_sample, led)

% Filters out short runs of data caused by tracker picking up an incorrect distant point.
% Resets led_pos_x to NaN

% Find first OK point & assume that this is OK
ok_pos = find( ~isnan(led_pos(:,led,1)) );

if length(ok_pos) < 2
    warning('scaNpix::dacqUtils::postprocess_posData::led_speed_filter: < 2 tracked points for led %d\n', led);
end
mpps_sqd = max_pix_per_sample^2;

n_jumpy  = 0;
prev_pos = ok_pos(1);
for i = 2:length(ok_pos)
    pos = ok_pos(i);
    % Get speed of shift from prev_pos in pixels per sample (squared)
    pix_per_sample_sqd = ((led_pos(pos,led,1)-led_pos(prev_pos,led,1))^2+(led_pos(pos,led,2)-led_pos(prev_pos,led,2))^2)/(pos-prev_pos)^2;
    if pix_per_sample_sqd > mpps_sqd
        led_pos(pos, led, 1) = NaN;
        n_jumpy              = n_jumpy+1;
    else
        prev_pos             = pos;
    end
end
end

% ----------------------------------------------------------------------------------------------------------------------

function swap_list = led_swap_filter(led_pos, led_pix)
% Checks for instances of two leds swapping or big one replacing little one
% when the big one gets obscured.
% Input xy posiiton of each led and
% and number of pixels in each. Big light is light number 1
% format: led_pos(1:n_pos, 1:num_cols, x-y), npix(1:n_pos, 1:num_cols)

thresh    = 5;

mean_npix = mean(led_pix,'omitnan');
std_npix  =  std(led_pix,[],'omitnan');

% Check if big light closer to small light at t-1 than to big light at t-1
% and small light either not found or closer to big light at t-1 than small light at t-1
pos       = 2:size(led_pix,1);

% Use one of the two following blocks of code - the first calculates a city
% block metric and the second euclidian distance. For most applications the
% latter is correct.

% Calculate city block metric
% dist12 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 2, 2));
% dist11 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 1, 2));
% dist21 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 1, 2));
% dist22 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 2, 2));

%Calculate eucldian
dist12   = sqrt(sum(((squeeze(led_pos(pos,1,:)) - squeeze(led_pos(pos-1,2,:))).^2)','omitnan'));
dist11   = sqrt(sum(((squeeze(led_pos(pos,1,:)) - squeeze(led_pos(pos-1,1,:))).^2)','omitnan'));
dist21   = sqrt(sum(((squeeze(led_pos(pos,2,:)) - squeeze(led_pos(pos-1,1,:))).^2)','omitnan'));
dist22   = sqrt(sum(((squeeze(led_pos(pos,2,:)) - squeeze(led_pos(pos-1,2,:))).^2)','omitnan'));

switched = (dist12 < dist11-thresh) & ( isnan(led_pos(pos, 2, 1))' | (dist21 < dist22-thresh) );

% Check if size of big light has shrunk to be closer to that of small light (as Z score)
z11       = (mean_npix(1) - led_pix(pos, 1))/std_npix(1);
z12       = (led_pix(pos, 1) - mean_npix(2))/std_npix(2);
shrunk    = z11 > z12;
swap_list = find( switched & shrunk' ) + 1;
end

