function [xy, dir, speed, dir_disp] = postprocess_pos_data_v2(posdata, max_speed, box_car, setfile_header, headPos )
% Perform postprocessing on position data
% See: scanpix.dacqUtils
%
% [xy, dir, speed] = scanpix.dacqUtils.postprocess_pos_data_v2(posdata, max_speed, box_car, setfile_header);
%
% posdata is in format of global tint{x}.pos, ie .led_pos(1:n_pos, 1:n_leds, 2),
%                                                .led_pix(1:n_pos, 1:n_leds)
%                                                .header
% ARGUMENTS
% a fifth argument headPos can be specified being a real number between 0 and 1
% inclusive. Indicating where the animal's head is located between the two LEDs (0
% being at the large LED, 1 being at the small LED and other values indicating
% positions between those. [default value = 0.5]
% 
%
% Returns
% xy - position (in pixels - ij format),
% dir - direction in degrees anticlock wisefrom x axis,
% speed - in cm/s (using pixels_per_metre and position_sample_rate in posfile header).
% times - in s, vector of time of each pos sample starts with e.g. [1/50, 2/50, ...] 
% jumpyPercent - %of pos points that were excluded because of too high speed
% n_leds - number of LEDs used during recording
% dir_disp - direction (degs anti-clock from x-axis) derived from heading. In case of
%           n_leds == 1 then dir_disp==di
%
% org. mTint function with a few small edits
%

%% This function could use some TLC and a bit more clarity, esp. given we are not using mTint

%%
n_pos = size(posdata.led_pos,1); % Don't use led_pix because this is absent from the older format
n_leds = size(posdata.led_pos,2);

pos_sample_rate = sscanf(scanpix.dacqUtils.getValue(posdata.header, 'sample_rate'),'%d');

% For 2 spot tracking, check for instances of swapping (often happens if one LED bright, one less bright).
if n_leds == 2 && not(isempty(posdata.led_pix)) % Only check if we actually have led_pix
    swap_list = led_swap_filter(posdata.led_pos, posdata.led_pix);
    dum = posdata.led_pos(swap_list, 1, :);
    posdata.led_pos(swap_list, 1, :) = posdata.led_pos(swap_list, 2, :);
    posdata.led_pos(swap_list, 2, :) = dum;
    dum = posdata.led_pix(swap_list, 1);
    posdata.led_pix(swap_list, 1) = posdata.led_pix(swap_list, 2);
    posdata.led_pix(swap_list, 2) = dum;
end

% Filter points for those moving impossibly fast and set led_pos_x to NaN
pix_per_metre = sscanf(scanpix.dacqUtils.getValue(posdata.header,'pixels_per_metre'),'%d');
max_pix_per_sample = max_speed*pix_per_metre/pos_sample_rate;
for led = 1:n_leds
    [n_jumpy, posdata.led_pos] = led_speed_filter( posdata.led_pos, max_pix_per_sample, led);
    if( n_jumpy > n_pos/3 )
        warning(sprintf(' %d/%d positions rejected for led %d\n', n_jumpy, n_pos, led));
    end
end

% Interpolate to replace all NaN led positions
% 1/12/09 AJ: I've made changes to the following lines to make the pos.xy
% output more robust. interp1 ignores NaNs at the endpoints of a vector. So
% two lines have been added to take care of these. missing_pos is formed to
% reject poses where the system has failed to record either x OR y
% coordinate and ok_pos requires both x AND y coords to be present.
for led = 1:n_leds
        missing_pos = find(isnan(posdata.led_pos(:, led, 1)) | isnan(posdata.led_pos(:, led, 2)));
        ok_pos = find(~isnan(posdata.led_pos(:, led, 1)) & ~isnan(posdata.led_pos(:, led, 2)));
    for k = 1:2
        posdata.led_pos(missing_pos, led, k) = interp1(ok_pos, posdata.led_pos(ok_pos, led, k), missing_pos, 'linear');
        posdata.led_pos(missing_pos(missing_pos>max(ok_pos)), led, k) = posdata.led_pos(max(ok_pos), led, k);
        posdata.led_pos(missing_pos(missing_pos<min(ok_pos)), led, k) = posdata.led_pos(min(ok_pos), led, k);
    end
end

% Estimate position, direction and speed from led_pos, using smoothing
if( box_car > 0 )
    b = ones(ceil(box_car*pos_sample_rate), 1)./ceil(box_car*pos_sample_rate);
else
    b = 1;
end

% Need to know angles (and distances) of LEDs from rat (in bearing_colour_i in .pos header)
pos = 1:n_pos;
pos2 = 1:n_pos-1;
xy = ones(n_pos, 2)*NaN;
if n_leds == 1
    xy(pos, :) = posdata.led_pos(pos, 1, :);
    xy = imfilter(xy, b, 'replicate');
    % y from tracker increases downwards. dir is positive anticlockwise from X axis, 0<= dir <360
    %CB verified calc of dir in next line is correct & concords with tint
    dir(pos2) = mod((180/pi)*(atan2( -xy(pos2+1, 2)+xy(pos2, 2), xy(pos2+1, 1)-xy(pos2, 1) )), 360);
    dir(n_pos) = dir(n_pos-1);
    dir_disp=dir; %Return dir_disp for completness even though == dir
elseif n_leds == 2
    % 2 LEDs are assumed to be at 180deg to each other with their midpoint over the
    % animals head. Convention is that colbearings_set(1) is the large light (normally
    % at front) and (2) is the small light. Position of lights (defined in set
    % file header) is defined in degs with 0 being straight ahead of animal, values
    % increase anti-clockwise.
    
    % Smooth lights individually, then get direction. %% TW, 12/09/08: This method replicates TINT.
    smLightFront(pos, :) = imfilter(posdata.led_pos(pos, 1, :), b, 'replicate');
    smLightBack(pos, :) = imfilter(posdata.led_pos(pos, 2, :), b, 'replicate');
    
    correction = setfile_header.lightBearing(1); %To correct for light pos relative to rat subtract angle of large light
    dir = mod((180/pi)*( atan2(-smLightFront(pos,2)+smLightBack(pos,2), +smLightFront(pos,1)-smLightBack(pos,1)) ) - correction, 360);
    dir = dir(:);
    % Get position from smoothed individual lights %%  % TW, 12/09/08
    wghtLightFront=1-headPos;
    wghtLightBack=headPos;
    xy(pos, :) = (smLightFront(pos, :).*wghtLightFront + smLightBack(pos, :).*wghtLightBack);  %%% CB added code todo headPos other than 0.5.
    clear wght* headPos
    
    %Get heading from displacement too
    dir_disp(pos2) = mod((180/pi)*(atan2( -xy(pos2+1, 2)+xy(pos2, 2), xy(pos2+1, 1)-xy(pos2, 1) )), 360);
    dir_disp(n_pos) = dir_disp(n_pos-1);
    dir_disp=dir_disp(:);
end

%%% Calculate speed, based on distance(sampleN+1-sampleN) %%%
speed(pos2) = sqrt((xy(pos2+1,1)-xy(pos2,1)).^2+(xy(pos2+1,2)-xy(pos2,2)).^2);
speed(n_pos) = speed(n_pos-1);
speed = speed.*(100*pos_sample_rate/pix_per_metre);
speed = speed(:);

% times = (1:n_pos)/pos_sample_rate;
% times = times(:);


function [n_jumpy, led_pos] = led_speed_filter(led_pos, max_pix_per_sample, led)

% Filters out short runs of data caused by tracker picking up an incorrect distant point.
% Resets led_pos_x to NaN

% Find first OK point & assume that this is OK
ok_pos = find( ~isnan(led_pos(:,led,1)) );

if length(ok_pos) < 2
    warning(sprintf(' < 2 tracked points for led %d\n', led));
end
mpps_sqd = max_pix_per_sample^2;

n_jumpy = 0;
prev_pos = ok_pos(1);
for i = 2:length(ok_pos)
    pos = ok_pos(i);
    % Get speed of shift from prev_pos in pixels per sample (squared)
    pix_per_sample_sqd = ((led_pos(pos,led,1)-led_pos(prev_pos,led,1))^2+...
        (led_pos(pos,led,2)-led_pos(prev_pos,led,2))^2)/(pos-prev_pos)^2;
    if pix_per_sample_sqd > mpps_sqd
        led_pos(pos, led, 1) = NaN;
        n_jumpy = n_jumpy+1;
    else
        prev_pos = pos;
    end
end

% ----------------------------------------------------------------------------------------------------------------------

function swap_list = led_swap_filter(led_pos, led_pix)
% Checks for instances of two leds swapping or big one replacing little one
% when the big one gets obscured.
% Input xy posiiton of each led and
% and number of pixels in each. Big light is light number 1
% format: led_pos(1:n_pos, 1:num_cols, x-y), npix(1:n_pos, 1:num_cols)

thresh = 5;

mean_npix = nanmean(led_pix);
std_npix = nanstd(led_pix);

% Check if big light closer to small light at t-1 than to big light at t-1
% and small light either not found or closer to big light at t-1 than small light at t-1
pos = 2:size(led_pix,1);

% Use one of the two following blocks of code - the first calculates a city
% block metric and the second euclidian distance. For most applications the
% latter is correct.

% Calculate city block metric
% dist12 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 2, 2));
% dist11 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 1, 2));
% dist21 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 1, 2));
% dist22 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 2, 2));

%Calculate eucldian
dist12 = sqrt(nansum(((squeeze(led_pos(pos,1,:))-squeeze(led_pos(pos-1,2,:))).^2)'));
dist11 = sqrt(nansum(((squeeze(led_pos(pos,1,:))-squeeze(led_pos(pos-1,1,:))).^2)'));
dist21 = sqrt(nansum(((squeeze(led_pos(pos,2,:))-squeeze(led_pos(pos-1,1,:))).^2)'));
dist22 = sqrt(nansum(((squeeze(led_pos(pos,2,:))-squeeze(led_pos(pos-1,2,:))).^2)'));

switched = (dist12 < dist11-thresh) & ( isnan(led_pos(pos, 2, 1))' |(dist21 < dist22-thresh) );

% Check if size of big light has shrunk to be closer to that of small light (as Z score)
z11 = (mean_npix(1) - led_pix(pos, 1))/std_npix(1);
z12 = (led_pix(pos, 1) - mean_npix(2))/std_npix(2);
shrunk = z11 > z12;
swap_list = find( switched & shrunk' ) + 1;

% ----------------------------------------------------------------------------------------------------------------------