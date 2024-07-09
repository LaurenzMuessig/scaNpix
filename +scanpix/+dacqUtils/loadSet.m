function loadSet(obj, trialIterator)
% loadSet - load set file data from DACQ files
%
% Syntax:  obj = loadSet(obj,trialIterator)
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

% sanity check
if ~exist(fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.set']),'file')
    ME = MException('scaNpix:loadSet:setFileNotFound', ['Could not find ''' strrep(obj.dataPath{trialIterator},'\','/') obj.trialNames{trialIterator} '.set''' ]); % should maybe switch to forward slash everywhere
    throw(ME);
end

% Read file %
fid = fopen(fullfile(obj.dataPath{trialIterator}, [obj.trialNames{trialIterator} '.set']), 'r');
C = textscan(fid, '%s %[^\r\n]');
sFileTxt = [cat(1,C{1}) cat(1,C{2})];
fclose(fid);

% legacy data isn't supported. multiplex pre-amp has different ADC range
if str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ADC_fullscale_mv')) ~= 1500 && isempty(scanpix.dacqUtils.getValue(sFileTxt, 'demux_en_dac_1'))
    ME = MException('scaNpix:loadSet:legacyDataError', ['''' strrep(obj.dataPath{trialIterator},'\','/') obj.trialNames{trialIterator} '.set'' was not acquired with DACQ USB. Data obtained from legacy versions are not supported.' ]);
    throw(ME);
end

%%% These are straightforward one value fields %%%
setFile.tracked_spots    = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'tracked_spots'));
setFile.xmin             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'xmin'));
setFile.xmax             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'xmax'));
setFile.ymin             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ymin'));
setFile.ymax             = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ymax'));
setFile.sw_version       = scanpix.dacqUtils.getValue(sFileTxt, 'sw_version');   % Don't STR2DOUBLE - can be 4.00a
setFile.trial_time       = scanpix.dacqUtils.getValue(sFileTxt, 'trial_time');
setFile.ADC_fullscale_mv = str2double(scanpix.dacqUtils.getValue(sFileTxt, 'ADC_fullscale_mv'));
duration                 = sscanf(scanpix.dacqUtils.getValue(sFileTxt, 'duration'),'%d');
if ~~rem(duration,10)
    duration             = duration - rem(duration,10);
end
setFile.duration         = duration;

% light parameters: make a 1:4 vector %
lightParams = {'lightBearing', 'colactive'};
for i = 1:2
    for j = 1:4
        setFile.(lightParams{i})(j) = str2double(scanpix.dacqUtils.getValue(sFileTxt, [lightParams{i} '_' num2str(j)]));
    end
end

% Gains %
for i = 1:128 % HARD CODED (but Jim will probably never make a higher channel count system...)
    setFile.gains(i) = str2double(scanpix.dacqUtils.getValue(sFileTxt, ['gain_ch_' num2str(i-1)]));
end
setFile.gains     = reshape(setFile.gains,4,32)';
setFile.fullscale = (setFile.ADC_fullscale_mv ./ setFile.gains) .* 1000;

%%% EEG Channels %%%
% Which channels EEGs recorded? %
recordingChannel = zeros([1 128]); %HARD CODED
for i = 1:length(recordingChannel) %%
    temp=scanpix.dacqUtils.getValue(sFileTxt,['saveEEG_ch_' num2str(i)]);
    if isempty(temp)
        break;
    elseif str2double(temp) % temp is '1' or '0' for EEG used or not.
        recordingChannel(i) = str2double(scanpix.dacqUtils.getValue(sFileTxt,['EEG_ch_'  num2str(i)]));
    end
end
EEGSlotActive    = find(recordingChannel);                  % Now a list of eeg channels in use ..
recordingChannel = recordingChannel(recordingChannel~=0);   %  .. and their recording slots. (recordingChannel=0 if slot not in use).

% Get signal source channel, gain, filters %
if isempty(recordingChannel)
    [setFile.lfp_channel,setFile.lfp_recordingChannel,setFile.lfp_slot,setFile.lfp_scalemax,setFile.lfp_filter,setFile.lfp_filtresp,setFile.lfp_filtkind,setFile.lfp_filtfreq1,...
        setFile.filtfreq2,setFile.lfp_filtripple] = deal([]); % In case of null EEG
else
    for i=1:length(recordingChannel)
        mode = str2double(scanpix.dacqUtils.getValue(sFileTxt, ['mode_ch_' num2str(recordingChannel(i)-1)]));
        if any( mode == [1 3] ) % Mode B or -B (ref or -ref)
            chTemp = scanpix.dacqUtils.getValue(sFileTxt, ['b_in_ch_' num2str(recordingChannel(i)-1)]);
            % if not free referencing system need to fetch ref channel ID from ref field
            if isempty(scanpix.dacqUtils.getValue(sFileTxt,'modeanalog32'))
                chTemp = scanpix.dacqUtils.getValue(sFileTxt, ['ref_' chTemp]);
            end
            setFile.lfp_channel(i) = str2double(chTemp) + 1;
        else
            setFile.lfp_channel(i) = recordingChannel(i);
        end
        setFile.lfp_recordingChannel(i) = recordingChannel(i);   % Actual recording channel, not signal source (useful when approaching raw data).
        setFile.lfp_slot(i)             = EEGSlotActive(i);      % Is this EEG1, EEG2, etc (according to the 'record>setup>EEG' tab.
        setFile.lfp_scalemax(i)         = setFile.fullscale( ceil(recordingChannel(i)/4), recordingChannel(i)-(((ceil(recordingChannel(i)/4))-1)*4) );
        setFile.lfp_filter(i)           = str2double(scanpix.dacqUtils.getValue(sFileTxt, ['filter_ch_' num2str(recordingChannel(i)-1)]));
        
        % DacqUSB only - fuller filter spec %
        fSpec = {'filtresp', 'filtkind', 'filtfreq1', 'filtfreq2', 'filtripple'};
        for j = 1:length(fSpec)
            setFile.(['lfp_' fSpec{j}]) = str2double(scanpix.dacqUtils.getValue(sFileTxt, [fSpec{j} '_ch_' num2str(recordingChannel(i)-1)]));
        end
        
    end
end
% need to initialise here
setFile.ppm             = []; % will be populated when loading pos data
setFile.ppm_org         = []; % will be populated when loading pos data
% for these the trick is to add them after initial loading using the addMetaData method and then reload the pos data with the reload flag
setFile.trialType       = []; % no method yet for dacq data to gather that
setFile.trackLength     = []; % no method yet for dacq data to gather that
setFile.envSize         = []; % no method yet for dacq data to gather that
setFile.envBorderCoords = []; % no method yet for dacq data to gather that

% output
if isempty(obj.trialMetaData)
    obj.trialMetaData = setFile;
else
    obj.trialMetaData(trialIterator) = setFile;
end

if isempty(obj.dataSetName)
    obj.dataSetName = ['r_' obj.trialNames{trialIterator}(1:6)]; %%%%% THIS NEEDS FIX ONCE WE ARE CLEAR HOW TO HANDLE METADATA IN DACQ
end

end

