function extract_waveforms(ephysObj,trialInd,varargin)
% extract_waveforms - extract waveform data from raw neuropixel data.
% Actual data is grabbed in subfunction.  
%
% Usage:
%       [npixObj,waveforms,channels] = scanpix.npixUtils.extract_waveforms(npixObj);
%       [npixObj,waveforms,channels] = scanpix.npixUtils.extract_waveforms(npixObj,trialInd);
%       [npixObj,waveforms,channels] = scanpix.npixUtils.extract_waveforms(__, Name-Value comma separated list);
%
%
% Inputs:   npixObj     - npix class object
%           trialInd    - index for trial we want to get waveform data for  
%           varargin    - name-value: comma separated list of name-value pairs     
%
% Outputs:  npixObj     - class obj, now containing waveform data
%           waveforms   - array of waveforms for each cluster (format:)
%           channels    - array of channel IDs for each cluster
%
% LM 2021
%
%% TO DO:

%% PARAMS
p = parseParams(varargin); % see subfunction for parameter list
% not sure this is still useful as we typically load waveforms for all
% trials
if nargin == 1 || isempty(trialInd) || strcmp(p.Results.mode,'cat')
    trialInd = 1:length(ephysObj.trialNames);
end

%
driftFlag = false;
nChan     = 385; %ephysObj.trialMetaData(1).nChanTot;

%% deal with raw data input
if strcmp(p.Results.mode,'cat')
    % load from drift corr file should be the default really whjen using
    % KS2.5 or 3
    if isempty(p.Results.path2cat)
        % prompt user to select file
        [fNameCat,pathCat] = uigetfile(fullfile(cd,'*.dat;*.ap.bin'),'Please Select the Drift Corrected File');
        if isnumeric(pathCat)
            warning('If you want to use the drift corrected data to extract waveforms, I need some info where that might be found on your disk. Too late now, but maybe you''ll do better later.');
            return
        end
        path2raw = {[pathCat fNameCat]};
    else
        path2raw = p.Results.path2cat;
        if ~iscell(path2raw)
            path2raw = {path2raw};
        end
    end
    % need a few details from concat log file
    fsLog   = dir(fullfile(fileparts(path2raw{1}),'*logFile.tsv'));
    logFile = tdfread([fsLog.folder filesep fsLog.name],'tab');
    %
    [~,~,ext] = fileparts(path2raw{1});
    if strcmp(ext,'.dat')
        driftFlag = true;
        nChan     = ephysObj.trialMetaData(1).nChanSort;
    end
elseif strcmp(p.Results.mode,'single')
    % if you load from ap.bin raw - you should really HP filter and CAR this data before extracting waveforms
    path2raw = fullfile(ephysObj.dataPath(trialInd),strcat(ephysObj.trialNames(trialInd),ephysObj.fileType));
else
    error([p.Results.mode ' is not a valid option. Try ''cat'' or ''single'' instead.']);
end

%% deal with clusters to extract
if ~isempty(p.Results.clu)
    cluInd = ismember( ephysObj.cell_ID(:,1),p.Results.clu);
else
    cluInd = true(length( ephysObj.cell_ID(:,1)),1);
end
spkCount = cumsum(cluInd);

%%
% if p.Results.peelTmpl
%     
% %     templates = readNPY(fullfile(fileparts(path2raw{1}),'templates.npy'));
% %     ephysObj.dataPathSort{1} = 'Z:\r1158_20230606T171103\r1158_220518_g1_t0\';
%     templates = readNPY(fullfile(fileparts(ephysObj.dataPathSort{1}),'templates.npy')); 
%     relTvec.temp = -floor( size(templates, 2)/2)+1:-floor( size(templates, 2)/2)+ size(templates, 2); % assumes template occurs at left edge
%     relTvec.snip = -floor(p.Results.nsamp/2):floor(p.Results.nsamp/2);
%     
% end

% ephysObj.dataPath{1} = 'S:\1postDoc\Neuropixels\rawData\r1158\220518\r1158_220518_g1\r1158_220518_g1_imec0\';
%% loop over trials and clusters
hWait = waitbar(0);
for i = trialInd % loop over trials
    % only read concat file once!
    if strcmp(p.Results.mode,'single') || (i == 1 && strcmp(p.Results.mode,'cat'))
        binFileStruct  = dir( path2raw{i} );
        if isempty(binFileStruct)
            error('Can''t find raw data mate...!')
        end
        dataTypeNBytes = numel(typecast(cast(0, p.Results.prec), 'uint8')); % determine number of bytes per sample
        nSamp          = binFileStruct.bytes/(nChan*dataTypeNBytes);  % Number of samples per channel - quicker than reading from meta file

        % we need to check in channel map if recording spanned multiple banks -
        % NEEDS WORK TO ACCOUNT BETTER FOR ALL EVENTUALITIES
        chanMapFile    = dir( fullfile(binFileStruct.folder, '*ChanMap.mat') );
        if isempty(chanMapFile)
            chanMapFName = scanpix.npixUtils.SGLXMetaToCoords_v2(binFileStruct.folder);
        else
            if driftFlag
                chanMapFile    = dir( fullfile(binFileStruct.folder, '*driftCorrChanMap.mat') );
                % in case there is no drift corr channel map (legacy data), we need
                % to generate it
                if isempty(chanMapFile)
                    tmpChanMap    = dir( fullfile(binFileStruct.folder, '*kilosortChanMap.mat') );
                    saveDriftCorrChanMap(fullfile(tmpChanMap.folder,tmpChanMap.name));
                    chanMapFile    = dir( fullfile(binFileStruct.folder, '*driftCorrChanMap.mat') );
                end
                chanMapFName = fullfile(chanMapFile.folder,chanMapFile.name);
            else
                chanMapFile    = dir( fullfile(binFileStruct.folder, '*kilosortChanMap.mat') );
                chanMapFName = fullfile(chanMapFile.folder,chanMapFile.name);
            end
        end
        chanMap        = load(chanMapFName);
        bankBoundaries = find(abs(diff(chanMap.ycoords)) > p.Results.chanspace) + [0 1];
        
        if driftFlag
            try
%                 tmp  = load(fullfile(fileparts(path2raw{i}),'whiteMat.mat')); % for KS3 we can't use the .npy anymore as it's just a diagonal matrix with the int16 scaling
                tmp  = load(fullfile(fileparts(ephysObj.dataPathSort{i}),'whiteMat.mat'));
                winv = tmp.Wrot^-1;
            catch
                error('Can''t find whitening matrix! Either find that file or load from raw .ap.bin. Does that sound like a plan?');
%                 winv = readNPY(fullfile(fileparts(path2raw{i}),'whitening_mat_inv.npy'));
%                 winv = 1;
            end
        else
            winv = 1;
        end
        
        nSampPrePeak = round(p.Results.prepeak*p.Results.nsamp);
        
%         if p.Results.peelTmpl
% %             templates = readNPY(fullfile(fileparts(path2raw{i}),'templates.npy'));
%             templates = readNPY(fullfile(fileparts(ephysObj.dataPathSort{i}),'templates.npy')); 
% %             winvTmp = readNPY(fullfile(fileparts(path2raw{i}),'whitening_mat_inv.npy'));
%             winvTmp = readNPY(fullfile(fileparts(ephysObj.dataPathSort{i}),'whitening_mat_inv.npy'));
%             sz = size(templates);
%             templatesUnWhite = reshape(reshape(templates, sz(1)*sz(2), sz(3)) * winvTmp, sz);
%             [~, resortInd] = sort(ephysObj.cell_ID(:,1)); 
%             templatesUnWhite = templatesUnWhite(resortInd,:,:);
%             %
% %             Z:\r1158_20230606T171103
% % %             amps = readNPY(fullfile(fileparts(path2raw{i}),'amplitudes.npy'));
% %             amps = readNPY(fullfile('Z:\r1158_20230606T171103\r1158_220518_g1_t0','amplitudes.npy'));
% % %             tmpST = readNPY(fullfile(fileparts(path2raw{i}),'spike_times.npy'));
% %             tmpST = readNPY(fullfile('Z:\r1158_20230606T171103\r1158_220518_g1_t0','spike_times.npy'));
%             %
%             STVect = vertcat(ephysObj.spikeData.spk_Times{i}{:});
%             nSPikesPerClu = cellfun(@length,ephysObj.spikeData.spk_Times{i});
%             cluVect = cell2mat(arrayfun(@(x,y) repmat(x,y,1),ephysObj.cell_ID(:,1),nSPikesPerClu,'uni',0));
%             
% %             ephysObj.dataPathSort{i} = 'Z:\r1158_20230606T171103\r1158_220518_g1_t0\';
%             amps = scanpix.npixUtils.loadSpikesNPix(ephysObj,i,0,1);
%             amps = vertcat(amps{:});
% %             templateIDs = vertcat(templateIDs{:});
%         end
        
        
        %% load file
        mmf = memmapfile(path2raw{i}, 'Format', {p.Results.prec, [nChan nSamp], 'x'});

    end
    
    %% extract waveforms    
    if driftFlag
        ind     = find(strcmp([ephysObj.trialNames{i} ephysObj.fileType],cellstr(logFile.filename)));
        tempST  = cellfun(@(x) x + ephysObj.trialMetaData(i).offSet + sum(logFile.duration(2:ind-1)), ephysObj.spikeData.spk_Times{i},'uni',0); % add offset to spike times
%         if p.Results.peelTmpl; STVect = ceil( (STVect + ephysObj.trialMetaData(i).offSet + sum(logFile.duration(2:ind-1))) .* ephysObj.params('APFs')); end
    else
        tempST  = cellfun(@(x) x + ephysObj.trialMetaData(i).offSet, ephysObj.spikeData.spk_Times{i},'uni',0); % add offset to spike times
%         if p.Results.peelTmpl; STVect = ceil( (STVect + ephysObj.trialMetaData(i).offSet) .*  ephysObj.params('APFs')); end
    end

    
    [tmpWaveforms, tmpChannels] = deal(cell(length(tempST),1));
   
    for j = 1:length(tempST) % loop over cells
        
        if ~cluInd(j);continue;end
        
        currSTimesBin = ceil(tempST{j} * p.Results.fs);
        currChannels  = max([1,ephysObj.cell_ID(j,3)-p.Results.getnch]):min([nChan, ephysObj.cell_ID(j,3)+p.Results.getnch]); % take care not to go <0 or > nCHan
        
        % remove channels from list
        ind = ismember(currChannels,p.Results.remchans);
        if any(ind)
            currChannels = currChannels(~ind);
            lhsAdd = sum(find(ind) <= p.Results.getnch);
            if lhsAdd > 0;  currChannels = [currChannels(1)-lhsAdd:currChannels(1)-1, currChannels]; end
            rhsAdd = sum(find(ind) > p.Results.getnch+1);
            if rhsAdd > 0;  currChannels = [currChannels, currChannels(end)+1:currChannels(end)+rhsAdd]; end
        end
        
        % need to remove channels if selection spans multiple banks
        if any(ismember(currChannels,bankBoundaries))
            if sum(currChannels <= bankBoundaries(1)) > sum(currChannels >= bankBoundaries(2))
                currChannels = currChannels(currChannels <= bankBoundaries(1));
            else
                currChannels = currChannels(currChannels >= bankBoundaries(2));
            end
        end
        
        % now extract waveforms for current cluster
        if ~isempty(p.Results.nwave)
            nWFs2Extract = min(length(currSTimesBin),p.Results.nwave);
        else
            nWFs2Extract = length(currSTimesBin);
        end
        
        if strcmp(p.Results.mode,'single')
            currChannels = scanpix.npixUtils.mapChans(chanMap.connected,currChannels);
        end
        
        currWave     = nan(nWFs2Extract,p.Results.nsamp+1,2*p.Results.getnch+1);
        ind2extract  = ceil(linspace(1,length(currSTimesBin),nWFs2Extract));
        c = 1;
        for k = ind2extract
            
            
            if p.Results.filter
                startIdx = max([currSTimesBin(k)-0.25*p.Results.fs,1]);
                endIdx   = min([currSTimesBin(k)+0.25*p.Results.fs,nSamp]);
                currData = double(mmf.Data.x(:,startIdx:endIdx)'); %* winv;
%                 currData = mmf.Data.x(:,startIdx:endIdx)'; %* winv;
                currData = HPfilter(currData, currChannels, p.Results.fs);
%                 currData = currData';
                currData = currData(0.25*p.Results.fs-nSampPrePeak+1:0.25*p.Results.fs+p.Results.nsamp-nSampPrePeak+1,:); 
            else
                startIdx = max([currSTimesBin(k)-nSampPrePeak,1]);
                endIdx   = min([currSTimesBin(k)+p.Results.nsamp-nSampPrePeak,nSamp]);
                currData = double(mmf.Data.x(:,startIdx:endIdx)') * winv;
%                 currData = currData(:,currChannels); 
            end
            %
%             currData = currData .* p.Results.bitres ./ p.Results.gain .* 1e6; %

%             if p.Results.peelTmpl
%                 times = [currSTimesBin(k) currSTimesBin(k) + relTvec.snip(1) - relTvec.temp(end) currSTimesBin(k) + relTvec.snip(end) - relTvec.temp(1)];
%                 otherSpikes = peelOffSpikes(times, templatesUnWhite, ephysObj.cell_ID(j,1), currChannels, STVect, cluVect, amps, relTvec);
%                 currData = double(int16(currData(:,currChannels)) - int16(otherSpikes'));
%             else
                currData = currData(:,currChannels);
%             end
            currWave(c,:,:) = currData .* p.Results.bitres ./ p.Results.gain .* 1e6; %
%             currWave(c,:,:) = currData(:,currChannels);
%             currData = currData(:,currChannels); 
            %
%             currWave(c,:,1:length(currChannels)) = currData .* p.Results.bitres ./ p.Results.gain .* 1e6; % uV conversion
            c = c+1;
        end
        
       tmpWaveforms{j} = currWave;
       tmpChannels{j}  = currChannels';
        
        waitbar( (spkCount(j)*i)/(sum(cluInd)*length(trialInd)), hWait, ['cluster ' num2str(spkCount(j)) '/' num2str(sum(cluInd)) ' from trial ' num2str(i) '/' num2str(length(trialInd))] );
    end
    %
    if isempty(ephysObj.spikeData.spk_waveforms{i})
        ephysObj.spikeData.spk_waveforms{i} = [tmpWaveforms, tmpChannels];
    else
        ind = ~cellfun('isempty',tmpWaveforms);
        ephysObj.spikeData.spk_waveforms{i}(ind,:) = [tmpWaveforms(ind),tmpChannels(ind)];
    end
    %
    if p.Results.save
        waveforms = [tmpWaveforms tmpChannels];
        save(fullfile(ephysObj.dataPath{i},'waveforms.mat'),'waveforms');
    end
    
end
close(hWait);

end

% subfunction to parse input
function p = parseParams(varargs)

defaultChanIgnore    = [];        % 192: reference channel; 385: sync channel - these should def be ignored
defaultPrecision     = 'int16';   % Data type of file
% defaultNCh           = 385;       % Number of channels in file
defaultNChWaves      = 5;         % grab +/- this many channels around peak channel of cluster
defaultNWaves        = 250;       % this many waveforms/cluster (if [] we'll grab all)
defaultNSamplesWf    = 69;        % this many samples/AP (40 = 1.3ms)
defaultPropSampPre   = 0.4;       % relative amount of samples pre peak (0.375 @ 40 samples = 15 samples pre peak and 34 samples post peak
defaultFs            = 30000;     % sampleing rate
defaultGain          = 500;       % gain
defaultBitResolution = 1.2/2^10;  % in V/bit
defaultChanSpaceVert = 20;        % vertical channel spacing in um
defaultApplyFilter   = false;     % do common average referencing
peelTemplates        = false;
% defaultUnWhitedata   = false;     % unwhiten data - important if using KS2.5 or higher as drift corrected file is whitened  
mode                 = 'single';
saveWFs              = false;
path2cat             = '';        % path to drift corr. file
cluIDs               = [];

p = inputParser;
addParameter(p,'remchans',   defaultChanIgnore);
addParameter(p,'prec',       defaultPrecision,@ischar);
% addParameter(p,'nch',        defaultNCh,@isscalar);
addParameter(p,'getnch',     defaultNChWaves,@isscalar);
addParameter(p,'nwave',      defaultNWaves,@(x) isscalar(x) || isempty(x));
addParameter(p,'nsamp',      defaultNSamplesWf,@isscalar);
addParameter(p,'prepeak',    defaultPropSampPre);
addParameter(p,'fs',         defaultFs,@isscalar);
addParameter(p,'gain',       defaultGain,@isscalar);
addParameter(p,'bitres',     defaultBitResolution,@isscalar);
addParameter(p,'chanspace',  defaultChanSpaceVert,@isscalar);
addParameter(p,'filter',     defaultApplyFilter);
addParameter(p,'peelTmpl',   peelTemplates, @islogical);
% addParameter(p,'unwhite',    defaultUnWhitedata);
addParameter(p,'save',       saveWFs,@islogical);
addParameter(p,'mode',       mode,@ischar);
addParameter(p,'path2cat',   path2cat,@(x) ischar(x) || iscell(x));
addParameter(p,'clu',        cluIDs);

parse(p,varargs{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveDriftCorrChanMap(pathChanMap)

tmp = load(pathChanMap);

chanMap     = (1:sum(tmp.connected))';
chanMap0ind = chanMap-1;
kcoords     = tmp.kcoords(tmp.connected);
xcoords     = tmp.xcoords(tmp.connected);
ycoords     = tmp.ycoords(tmp.connected);
connected   = true(sum(tmp.connected),1);

p = fileparts(pathChanMap);
[~,fn,~] = fileparts(tmp.name);
name = fullfile(p,[fn '_driftCorrChanMap.mat']);

save( name, 'chanMap', 'chanMap0ind', 'connected', 'name', 'xcoords', 'ycoords', 'kcoords' );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filtData = HPfilter(data, channels, fs)
% this should replicate the Kilosort HP filter (basically the fucnction
% 'gpuFilter' from the kilosort repo

fsHigh = 300;
fsLow  = 8000;

% set up the parameters of the filter
[b1, a1] = butter(3, [fsHigh/fs,fsLow/fs]*2, 'bandpass'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)

% subtract the mean from each channel
data = data - mean(data, 1); % subtract mean of each channel

% CAR, common average referencing by median
data = data - median(data, 2); % subtract median across channels

% next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)
data(:,channels) = filter(b1, a1, data(:,channels)); % causal forward filter
data(:,channels) = flipud(data(:,channels)); % reverse time
data(:,channels) = filter(b1, a1, data(:,channels)); % causal forward filter again
filtData = flipud(data); % reverse time back

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function otherSpikeData = peelOffSpikes(tWin, templates, currClu, currChans, spike_times, clusterIDs, amps, relTvec)
% % this should replicate the Kilosort HP filter (basically the fucnction
% % 'gpuFilter' from the kilosort repo
% 
% 
% %
% spikeInd = find(spike_times >= tWin(2) & spike_times <= tWin(3));  
% spikeInd(ismember(clusterIDs(spikeInd), currClu)) = [];
% %
% otherSpikeData = zeros(length(currChans),length(relTvec.snip));
% % loop over the nearby sp
% for i = 1:numel(spikeInd)
%     ind = spikeInd(i);
%     
%     % figure out time overlap and add to reconstruction
%     tOth = spike_times(ind);
%     indFromTemplate = relTvec.temp + tOth >= tWin(1) + relTvec.snip(1) & relTvec.temp + tOth <= tWin(1) + relTvec.snip(end);
%     indInsert = relTvec.snip + tWin(1) >= relTvec.temp(1) + tOth & relTvec.snip + tWin(1) <= relTvec.temp(end) + tOth;
%     insert = amps(ind) .* permute(templates(clusterIDs(ind), indFromTemplate, currChans), [3 2 1]);
%     if any(insert(:) > 0)
%         t = 1;
%     end
%     otherSpikeData(:,indInsert) = otherSpikeData(:,indInsert) + insert; 
% end
% 
% end



%                 
%               
% find spikes that would lie within this window (with padding), 
%                 % excluding those from clusters we wish to exclude
%                 t = int64(times(iT));
%                 minT = t + window(1) - int64(relTvec_template(end));
%                 maxT = t + window(2) - int64(relTvec_template(1));
%                 if iscell(exclude_cluster_ids_each_snippet)
%                     exclude_this = exclude_cluster_ids_each_snippet{iT};
%                 elseif ~isempty(exclude_cluster_ids_each_snippet)
%                     exclude_this = exclude_cluster_ids_each_snippet(iT);
%                 else
%                     exclude_this = [];
%                 end
%                 nearby_spike_inds = find(ds.spike_times >= minT & ds.spike_times <= maxT);      
%                 %nearby_spike_inds = find(ds.spike_times == t);
%                 
%                 nearby_spike_inds(ismember(ds.spike_clusters(nearby_spike_inds), exclude_this)) = [];
% 
%                 % figure out what channels we need to reconstruct onto
%                 cluster_ids_this = cluster_ids(iT);
%                 [~, cluster_ind_this] = ismember(cluster_ids_this, unique_cluster_ids);
%                 assert(cluster_ind_this > 0, 'Cluster for times(%d) not found in unique_cluster_ids', iT);
%                 channels_idx_this = channel_ids_by_cluster(:, cluster_ind_this);
%                 [~, channel_ind_this] = ismember(channels_idx_this, ds.channel_ids);
%                 assert(all(channel_ind_this) > 0, 'Some channels in channel_ids_by_cluster not found in channel_ids');
%                 
%                 if showPlots
%                     clf;
%                     plot(relTvec_snippet, p.Results.rawData(1, :, iT), 'k-', 'LineWidth', 2);
%                     hold on;
%                 end
%                 
%                 % loop over the enarby sp
%                 for iS = 1:numel(nearby_spike_inds)
%                     ind = nearby_spike_inds(iS);
%                     amp = ds.amplitudes(ind);
%                     
%                     % figure out time overlap and add to reconstruction
%                     tprime = int64(ds.spike_times(ind));
%                     indFromTemplate = relTvec_template + tprime >= t + relTvec_snippet(1) & relTvec_template + tprime <= t + relTvec_snippet(end);
%                     indInsert = relTvec_snippet + t >= relTvec_template(1) + tprime & relTvec_snippet + t <= relTvec_template(end) + tprime;
%                     insert = amp .* permute(templates(ds.spike_templates(ind), indFromTemplate, channel_ind_this), [3 2 1]);
%                     reconstruction(:, indInsert, iT) = reconstruction(:, indInsert, iT) + insert;
%                     
%                     if showPlots
%                         if tprime == t
%                             args = {'LineWidth', 1, 'Color', 'g'};
%                         else
%                             args = {};
%                         end
%                         plot(relTvec_template + tprime-t, amp .* templates(ds.spike_templates(ind), :, channel_ind_this(1)), 'Color', [0.7 0.7 0.7], args{:});
%                     end
%                 end
%                 
%                 if showPlots
%                     plot(relTvec_snippet, reconstruction(1, :, iT), 'r--', 'LineWidth', 2);
%                     plot(relTvec_snippet, p.Results.rawData(1, :, iT) - int16(reconstruction(1, :, iT)), 'b--');
%                     pause;
%                 end
%             end






