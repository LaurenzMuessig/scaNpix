function cell_ID = loadSpikes(obj, trialIterator)
% loadSpikes - load tetrode data from DACQ files
% We will load spike times and waveforms for each cluster
%
% Syntax:  cell_ID = loadSpikes(obj, trialIterator)
%
% Inputs:
%    obj           - ephys class object ('dacq')
%    trialIterator - numeric index for trial to be loaded
%
% Outputs:
%    cell_ID - numeric array of cell and tetrode IDs
%
% TW/LM 2020 (adapted from org. SCAN function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To-Do
% load MUA when no cut file is found
%
%%
[~,nTetrodes] = scanpix.dacqUtils.findDACQFiles( [obj.dataPath{trialIterator} filesep], obj.trialNames{trialIterator}, 'tet' );
if isempty(nTetrodes)
    ME = MException('scaNpix::loadSpikes:TetrodeFilesMissing', [' Can''t find any tetrode files for ' fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.set']) '.' ]);
    throw(ME);
end

clear scanpix.fxchange.textprogressbar  % get rid of persistent variable in 'textprogressbar' which can cause issues if you interrupt code during execution

% initialise
cell_ID                 = zeros(0, 2);
[spikeTimes, waveforms] = deal({});
scanpix.fxchange.textprogressbar(['Loading tetrode data for ' obj.trialNames{trialIterator} ' ']);

% loop through recorded tetrodes
for i = 1:length(nTetrodes)
    
    %%%%%%%%%%%%%% TETRODE FILE %%%%%%%%%%%%
    fid          = fopen(fullfile(obj.dataPath{trialIterator},[obj.trialNames{trialIterator} '.' num2str(nTetrodes(i))]),'r','ieee-be');        % Open tet file. 'ieee-be' string is machine format, which needs to be 'big endian' to read times correctly. Default format doesnt work on my PC.
    data         = fread(fid,'int8'); % Read data (header and voltage samples). 'int8' converts each byte into a integer -126:126.
    
    % extract some info from header
    hdr          = char(abs(data(1:400)))';                % Translate header into text. Use 'abs' to avoid error message on negative integers after data starts.
    % n spikes
    nspk_ind     = strfind(hdr,'num_spikes ') + 11;        % Look for number of spikes marker
    nspk         = sscanf(hdr(nspk_ind:end), '%d');
    % Get sample rate (DSP timer Ticks/s)
    timebase_ind = strfind(hdr,'timebase ') + 9;
    timebase     = sscanf(hdr(timebase_ind:end), '%d');
    
    % read actual data
    ds           = strfind(hdr,'data_start') + 10;         % Look for data start marker
    % spike data
    spk_AmplRaw  = data( ds : ds + nspk*216 - 1);          % Take voltage sample values. (Time stamps still included)
    spk_AmplRaw  = reshape(spk_AmplRaw,54,4,nspk);         % Reshape into (sample,channel,spike)
    spk_AmplRaw  = spk_AmplRaw(5:end,:,:);                 % Cut off time stamp samples
    %%%%%%
    
    % time stamps
    fseek(fid,(ds-1),'bof');                               % Set file position indicator at data start (for reading time stamps, see below)
    spk_TimesRaw = fread(fid,nspk,'int32',212);            % read time stamps. 'nspk'=read this many samples, 'int32'=read 32bit integer, '212'=skip 212 bytes between reading values.
    spk_TimesRaw = spk_TimesRaw ./ timebase;
    %%%%%%
    fclose(fid);
    
    %%%%%%%%%%%%%% CUT FILE %%%%%%%%%%%%
    % read cut file
    if strcmp(obj.params('cutFileType'),'cut')
        cutFileName = fullfile(obj.dataPath{trialIterator},[ obj.trialNames{trialIterator} obj.params('cutTag1') '_' num2str(nTetrodes(i)) obj.params('cutTag2') '.cut' ] );
    elseif strcmp(obj.params('cutFileType'),'clu')
        cutFileName = fullfile(obj.dataPath{trialIterator}, [ obj.trialNames{trialIterator} obj.params('cutTag1') '.clu.' num2str(nTetrodes(i)) ] );
    else
        ME = MException('scaNpix::loadSpikes:InvalidCutFileExtension', ['''' obj.params('cutFileType') '''' ' is not a valid extension for cut files - only ''cut'' or ''clu'' are allowed' ]);
        throw(ME);
    end
    
    if ~exist(cutFileName,'file')
        warning(['scaNpix: Can''t find ''' obj.trialNames{trialIterator} obj.params('cutTag1') '_' num2str(nTetrodes(i)) obj.params('cutTag2') '.cut''. No data loaded for Tet #' num2str(nTetrodes(i)) '.']);
        continue;
    end
    
    fid = fopen(cutFileName,'r');
    
    if strcmp(obj.params('cutFileType'),'cut')
        % get line with end of header (can be variable)
        lineN = 0;
        while 1
            tline = fgetl(fid);                                 % Look for start of exact cut data.
            if strncmp('Exact_cut_for',tline,12), lineN = lineN+1; break, end      %
            lineN = lineN+1;                                              % The final returned tline is the header for this.
        end
    else
        lineN = 1; % n.b. first line in clu files is n of clusters found
    end
    %
    frewind(fid);
    cut = cell2mat( textscan(fid,'%n','headerlines',lineN,'delimiter','/b') );
    fclose(fid);
    
    % sanity check
    if length(spk_TimesRaw) ~= length(cut)
        warning(['scaNpix: Trial ', obj.trialNames{trialIterator}, ': Number of spikes in cut and tetrode file for tetrode ' num2str(nTetrodes(i)) ' don''t match. Cut file corrupt?. Data loading for this tetrode will be skipped']);
        continue;
    end
    
    % remove spikes in overhang
    overhangInd  = spk_TimesRaw < obj.trialMetaData(trialIterator).duration;
    cut          = cut(overhangInd); % remove spikes in overhang
    spk_TimesRaw = spk_TimesRaw(overhangInd);
    spk_AmplRaw  = spk_AmplRaw(:,:,overhangInd);
    
    clustIds = unique(cut(cut~=0));  %  Note that .clu files have no convention of a cluster 0 that contains unsorted spikes (i.e. MUA and noise) - maybe we should make one?
    
    if isempty(clustIds)
        scanpix.fxchange.textprogressbar(i/length(nTetrodes)*100);
        continue;
    end
    
    [cut, sortInd]           = sort(cut);               % this will ensure output of accumarray is sorted properly (spike times in each cluster)
    spk_TimesRaw             = spk_TimesRaw(sortInd);   % sort other arrays as well accordingly
    spk_AmplRaw              = spk_AmplRaw(:,:,sortInd);
    % get spike times for each cluster
    temp_ST                  = accumarray(cut(cut~=0), spk_TimesRaw(cut~=0,1), [max(unique(cut)) 1], @(x) {x});
    ind_empty                = ~cellfun('isempty',temp_ST);
    temp_ST                  = temp_ST(ind_empty);
    % bit more complicated for waveform amplitudes
    % first convert to uV...
    chGains                  = repmat(obj.trialMetaData(trialIterator).fullscale(nTetrodes(i),:),[50,1,sum(overhangInd)]);
    spk_AmplRaw1             = spk_AmplRaw ./128 .* chGains;
    % ... then assign to clusters
    temp                     = shiftdim(repmat(cut, 1, 50, 4), 1); % all cluster IDs for raw amplitude array
    temp                     = temp(:);
    temp2                    = spk_AmplRaw1(:);
    temp_SA                  = accumarray(temp(temp~=0), temp2(temp~=0), [max(unique(temp)) 1],@(x) {x});
    temp_SA                  = cellfun(@(x) shiftdim( reshape(x,50,4,[] ), 2),temp_SA, 'UniformOutput',0); % restore old format ('reshape') & bring into same format as 'temp_ST' ('shiftdim')
    temp_SA                  = temp_SA(ind_empty);
    % accumulate output data
    spikeTimes(end+1:end+length(temp_ST),1) = temp_ST;
    waveforms(end+1:end+length(temp_ST),1)  = temp_SA;
    cell_ID(end+1:end+length(temp_ST),:)    = [clustIds, nTetrodes(i) * ones(length(clustIds),1)]; % [cell n, tet n]
    
    scanpix.fxchange.textprogressbar(i/length(nTetrodes)*100);
end
obj.spikeData(1).spk_Times{trialIterator}     = spikeTimes;
obj.spikeData(1).spk_waveforms{trialIterator} = waveforms;
obj.spikeData(1).sampleT{trialIterator}       = [];  % not relevant for dacq data, but we need to keep format consistent with neuropix data

scanpix.fxchange.textprogressbar('  DONE!');
end

