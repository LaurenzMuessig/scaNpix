function loadLFPs(obj, trialIterator)
% loadLFPs - load eeg data from DACQ files
% We will load all low and (if available) high sample rate EEGs and we also fetch the tetrode ID each eeg was recorded from
%
% Syntax:  loadLFPs(obj, trialIterator)
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
clear scanpix.fxchange.textprogressbar  % get rid of persistent variable in 'textprogressbar' which can cause issues if you stop code during execution

noEGFflag = false;
lfp2load = scanpix.dacqUtils.findDACQFiles([obj.dataPath{trialIterator} filesep], obj.trialNames{trialIterator},'eeg');
if ischar(lfp2load)
    lfp2load = {lfp2load}; % in case only 1 eeg
end
% sometimes useful to not load the high Fs eeg to save time
if obj.params('loadHighFsLFP')
    it = 1:2;
else
    it = 1;
end

sRateStr = {'lfpFs','lfpHighFs'};

for i = it
    % switch extension
    if i == 2
        lfp2load = strrep(lfp2load,'.eeg','.egf');
    end
    
    scanpix.fxchange.textprogressbar(['Loading ' lfp2load{1}(end-3:end) ' data for ' obj.trialNames{trialIterator} ' ']);
    
    for j = 1:length(lfp2load)
        % Note in some rare cases it seems that fopen gets the encoding scheme wrong and reads in a gibberish header so we want to be explicit 
        fid = fopen(fullfile(obj.dataPath{trialIterator}, lfp2load{j}),'r','ieee-be',"UTF-8");  % 'ieee-be' is machine format, 'big endian'.
        if fid == -1
            noEGFflag = true;
            scanpix.fxchange.textprogressbar(0);
            continue
        end
        
        hdr = fread(fid,400,'int8');
        ds  = strfind(hdr','data_start') + 10; % data start marker
        % more convenient format - read some info from header
        frewind(fid);
        %%%%%%%
        tempHeader = textscan(fid,'%s %[^\r\n]',11);
        tempHeader = horzcat(tempHeader{:});
        
        %%%%%%%
        sRateInd = strcmp('sample_rate',tempHeader(:,1));
        if sum( sRateInd ) == 0   % Sometimes there are 'empty' EEG files, which don't even have a full header. If we have one of these, just quit now.
            warning(['scaNpix: Problem loading EEG. Empty EEG data file for ' lfp2load{j}]);
            continue
        end
        obj.params(sRateStr{i}) = sscanf(tempHeader{sRateInd,2},'%d');  % hard code??
        %%%%%%
        bytesInd     = strcmp('bytes_per_sample',tempHeader(:,1));
        bytesPerSamp = sscanf(tempHeader{bytesInd,2},'%d');
        % read data
        fseek(fid,ds,'bof');
        if bytesPerSamp == 1
            nSamplesInd    = strcmp('num_EEG_samples',tempHeader(:,1));
            nSamples       = sscanf(tempHeader{nSamplesInd,2},'%d');
            tempData       = fread(fid,nSamples,'int8');
            obj.lfpData(1).lfp{trialIterator}{j}  = (double(tempData)./2^7) .* obj.trialMetaData(trialIterator).lfp_scalemax(j); %voltages
        elseif bytesPerSamp == 2
            nSamplesInd    = strcmp('num_EGF_samples',tempHeader(:,1));
            nSamples       = sscanf(tempHeader{nSamplesInd,2},'%d');
            %grab actual voltage data
            tempData = fread(fid,nSamples,'int16'); %re-read as int16
            obj.lfpData(1).lfpHighSamp{trialIterator}{j}  = (double(tempData)./2^15) .* obj.trialMetaData(trialIterator).lfp_scalemax(j); %voltages
        end
        fclose(fid);
        
        scanpix.fxchange.textprogressbar(j/length(lfp2load)*100);
    end
    if noEGFflag
        scanpix.fxchange.textprogressbar('  NO DATA!');
    else
        scanpix.fxchange.textprogressbar('  DONE!');
    end
end

obj.lfpData(1).lfpTet{trialIterator} = ceil( obj.trialMetaData(trialIterator).lfp_channel(:) ./ 4)';
end

