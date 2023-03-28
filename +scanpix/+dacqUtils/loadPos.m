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
nLights      = sum(obj.trialMetaData(trialIterator).colactive);
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
    posHeader{strcmp(posHeader(:,1),'pixels_per_metre'),2} = num2str(obj.params('ScalePos2PPM'));
    obj.trialMetaData(trialIterator).xmin = obj.trialMetaData(trialIterator).xmin * scaleFact;
    obj.trialMetaData(trialIterator).xmax = obj.trialMetaData(trialIterator).xmax * scaleFact;
    obj.trialMetaData(trialIterator).ymin = obj.trialMetaData(trialIterator).ymin * scaleFact;
    obj.trialMetaData(trialIterator).ymax = obj.trialMetaData(trialIterator).ymax * scaleFact;
end
obj.trialMetaData(trialIterator).ppm = sscanf(scanpix.dacqUtils.getValue(posHeader,'pixels_per_metre'),'%d');

% post process
posFile.header  = posHeader; % For mTint function compatibility.
posFile.led_pos = led_pos;
posFile.led_pix = led_pix;
%Interpolate, Filter and Smooth %  %%%%%%%%%% THIS COULD USE A BIT MORE EFFICIENCY AND BE UPDATED
[obj.posData(1).XYraw{trialIterator}, obj.posData(1).direction{trialIterator}, obj.posData(1).speed{trialIterator}] = scanpix.dacqUtils.postprocess_pos_data_v2(posFile, obj.params('posMaxSpeed'), obj.params('posSmooth'), obj.trialMetaData(trialIterator), obj.params('posHead'));

obj.params('posFs') = sscanf(scanpix.dacqUtils.getValue(posHeader,'sample_rate'),'%d');
% remove DACQ overhang
if obj.trialMetaData(trialIterator).duration * obj.params('posFs') < length(led_pos)
    obj.posData(1).XYraw{trialIterator}     = obj.posData(1).XYraw{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs'),:); % truncate data
    obj.posData(1).direction{trialIterator} = obj.posData(1).direction{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs')); % truncate data
    obj.posData(1).speed{trialIterator}     = obj.posData(1).speed{trialIterator}(1:obj.trialMetaData(trialIterator).duration * obj.params('posFs')); % truncate data
end

% convert to integers
obj.posData(1).XY{trialIterator} = [double( floor(obj.posData(1).XYraw{trialIterator}(:,1)) + 1 ), double( floor(obj.posData(1).XYraw{trialIterator}(:,2)) + 1 )];  %%% NECESSARY??

fprintf('  DONE!\n');

end

