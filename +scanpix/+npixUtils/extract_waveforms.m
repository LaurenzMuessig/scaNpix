function extract_waveforms(ephysObj,trialInd,options)
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
%           options     - name-value: comma separated list of name-value pairs (see arguments block)
%
% Outputs:
%
% LM 2021
%
%% TO DO:

%% PARAMS
arguments
    ephysObj
    trialInd                               = [];
    options.remchans  {mustBeNumeric}      = [];         % 192: reference channel; 385: sync channel - these should def be ignored
    options.prec      (1,:) char           = 'int16';    % Data type of file
    options.getnch    (1,1) {mustBeNumeric} = 5;         % grab +/- this many channels around peak channel of cluster
    options.nwave     {mustBeScalarOrEmpty} = 250;       % this many waveforms/cluster (if [] we'll grab all)
    options.nsamp     (1,1) {mustBeNumeric} = 69;        % this many samples/AP (40 = 1.3ms)
    options.prepeak   (1,1) {mustBeNumeric} = 0.4;       % relative amount of samples pre peak (0.375 @ 40 samples = 15 samples pre peak and 34 samples post peak
    options.fs        (1,1) {mustBeNumeric} = 30000;     % sampleing rate
    options.gain      (1,1) {mustBeNumeric} = 500;       % gain
    options.bitres    (1,1) {mustBeNumeric} = 1.2/2^10;  % in V/bit
    options.chanspace (1,1) {mustBeNumeric} = 20;        % vertical channel spacing in um
    options.filter    (1,1) logical         = false;     % do common average referencing and HP filtering
    options.fhigh     (1,1) {mustBeNumeric} = 300;       % filter high pass cut-off (Hz) - default matches CatGT -apfilter=butter,12,300,9000
    options.flow      {mustBeScalarOrEmpty} = 9000;      % filter low pass cut-off (Hz); [] for high pass only (kilosort default)
    options.forder    (1,1) {mustBeNumeric} = 3;         % butterworth order (applied forward + backward, as in kilosort)
    options.peelTmpl  (1,1) logical         = false;     % subtract templates of other clusters' spikes from snippets (drift corr. file only)
    options.tmplAlign {mustBeScalarOrEmpty} = [];        % template sample that corresponds to spike time; [] = from templates (median trough position)
    options.save      (1,1) logical         = false;
    options.mode      (1,:) char {mustBeMember(options.mode,{'single','cat'})} = 'single';
    options.path2cat  {mustBeA(options.path2cat,{'char','string','cell'})} = ''; % path to drift corr. file
    options.clu       {mustBeNumeric}      = [];
end

% not sure this is still useful as we typically load waveforms for all
% trials
if isempty(trialInd) || strcmp(options.mode,'cat')
    trialInd = 1:length(ephysObj.trialNames);
end

%
driftFlag = false;

%% deal with raw data input
if strcmp(options.mode,'cat')
    % load from drift corr file should be the default really when using KS2.5 or 3
    if isempty(options.path2cat)
        % prompt user to select file
        [fNameCat,pathCat] = uigetfile(fullfile(cd,'*.dat;*.ap.bin'),'Please Select the Drift Corrected File');
        if isnumeric(pathCat)
            warning('scaNpix::ephys::extract_waveforms:If you want to use the drift corrected data to extract waveforms, I need some info where that might be found on your disk. Too late now, but maybe you''ll do better later.');
            return
        end
        path2raw = {[pathCat fNameCat]};
    else
        path2raw = options.path2cat;
        if ~iscell(path2raw)
            path2raw = {path2raw};
        end
    end
    % need a few details from concat log file
    fsLog   = dir(fullfile(fileparts(path2raw{1}),'*logFile.tsv'));
    logFile = tdfread(fullfile(fsLog.folder,fsLog.name),'tab');
    %
    [~,~,ext] = fileparts(path2raw{1});
    if strcmp(ext,'.dat')
        driftFlag = true;
        nChan     = ephysObj.trialMetaData(1).nChanSort;
    else
        nChan     = ephysObj.trialMetaData(1).nChanTot;
    end
else % 'single' (mode is validated in arguments block)
    % if you load from ap.bin raw - you should really HP filter and CAR this data before extracting waveforms
    nChan    = ephysObj.trialMetaData(1).nChanTot;
    path2raw = fullfile(ephysObj.dataPath(trialInd),strcat(ephysObj.trialNames(trialInd),ephysObj.fileType));
end

% drift corr. file is already HP filtered and CAR'd (and whitened) by kilosort - filtering again would also skip unwhitening
applyFilter = options.filter;
if driftFlag && applyFilter
    warning('scaNpix::ephys::extract_waveforms:Drift corrected data is already filtered by kilosort. Ignoring ''filter'' = true.');
    applyFilter = false;
end

%% deal with clusters to extract
if ~isempty(options.clu)
    cluInd = ismember( ephysObj.cell_ID(:,1),options.clu);
else
    cluInd = true(length( ephysObj.cell_ID(:,1)),1);
end
spkCount   = cumsum(cluInd);

%% template peeling (subtract other clusters' spikes, as in neuropixel-utils 'subtractOtherClusters')
% only for drift corr. file - KS templates live in the same (filtered, whitened) space as that data
peelFlag = options.peelTmpl;
if peelFlag && ~driftFlag
    warning('scaNpix::ephys::extract_waveforms:Template peeling only works when extracting from the drift corrected .dat file. Ignoring ''peelTmpl'' = true.');
    peelFlag = false;
end

%% loop over trials and clusters
hWait = waitbar(0);
for i = trialInd % loop over trials
    % only read concat file once!
    if strcmp(options.mode,'single') || (i == 1 && strcmp(options.mode,'cat'))
        binFileStruct  = dir( path2raw{i} );
        if isempty(binFileStruct)
            error('scaNpix::ephys::extract_waveforms:Can''t find raw data mate...!')
        end
        dataTypeNBytes = numel(typecast(cast(0, options.prec), 'uint8')); % determine number of bytes per sample
        nSamp          = binFileStruct.bytes/(nChan*dataTypeNBytes);  % Number of samples per channel - quicker than reading from meta file

        % we need to check in channel map if recording spanned multiple banks -
        % NEEDS WORK TO ACCOUNT BETTER FOR ALL EVENTUALITIES
        chanMapFile = dir( fullfile(binFileStruct.folder, '*ChanMap.mat') );
        if isempty(chanMapFile)
            chanMapFName = scanpix.npixUtils.SGLXMetaToCoords_v2(binFileStruct.folder);
        else
            if driftFlag
                chanMapFile = dir( fullfile(binFileStruct.folder, '*driftCorrChanMap.mat') );
                % in case there is no drift corr channel map (legacy data), we need
                % to generate it
                if isempty(chanMapFile)
                    tmpChanMap  = dir( fullfile(binFileStruct.folder, '*kilosortChanMap.mat') );
                    saveDriftCorrChanMap(fullfile(tmpChanMap.folder,tmpChanMap.name));
                    chanMapFile = dir( fullfile(binFileStruct.folder, '*driftCorrChanMap.mat') );
                end
                chanMapFName = fullfile(chanMapFile.folder,chanMapFile.name);
            else
                chanMapFile  = dir( fullfile(binFileStruct.folder, '*kilosortChanMap.mat') );
                chanMapFName = fullfile(chanMapFile.folder,chanMapFile.name);
            end
        end
        chanMap        = load(chanMapFName);
        % work in sorted channel space (rows of KS/drift corr. data), same as ephysObj.cell_ID(:,3)
        sortedYCoords  = chanMap.ycoords(logical(chanMap.connected));
        nChanSorted    = numel(sortedYCoords);
        bankBoundaries = find(abs(diff(sortedYCoords(:))) > options.chanspace) + [0 1];
        
        if driftFlag
            try
                tmp  = load(fullfile(ephysObj.dataPathSort{i},'whiteMat.mat'));
                winv = tmp.Wrot^-1;
            catch
                error('scaNpix::ephys::extract_waveforms:Can''t find whitening matrix! Either find that file or load from raw .ap.bin. Does that sound like a plan?');
            end
        else
            winv = 1;
        end
        
        nSampPrePeak = round(options.prepeak*options.nsamp);
        wRel         = -nSampPrePeak:options.nsamp-nSampPrePeak; % snippet samples rel. to spike time

        %% load file
        mmf = memmapfile(path2raw{i}, 'Format', {options.prec, [nChan nSamp], 'x'});

        % sorting output of concatenated data lives with the drift corr. file (spike times are in cat file samples)
        if peelFlag
            tmpl = loadTemplates(fileparts(path2raw{i}), mmf, winv, nSamp, options.tmplAlign);
        end

    end
    
    %% extract waveforms    
    % spike times were split from the concatenated sort in samples, so trial offset in the cat file needs to be in exact samples
    % too (logFile.duration is fileTimeSecs = nSamp/imSampRate, which doesn't match nSamp/APFs). Row 1 of log is the cat file itself
    if strcmp(options.mode,'cat')
        ind        = find(strcmp([ephysObj.trialNames{i} ephysObj.fileType],cellstr(logFile.filename)));
        if isempty(ind)
            error(['scaNpix::ephys::extract_waveforms:Can''t find ' ephysObj.trialNames{i} ' in concat log file.']);
        end
        sampOffset = sum(logFile.nSamples(2:ind-1) ./ logFile.nChan(2:ind-1));
    else
        sampOffset = 0;
    end
    tempST  = cellfun(@(x) x + ephysObj.trialMetaData(i).offSet, ephysObj.spikeData.spk_Times{i},'uni',0); % add offset to spike times (in s, rel. to trial file start)

    %
    [tmpWaveforms, tmpChannels] = deal(cell(length(tempST),1));
   
    for j = 1:length(tempST) % loop over cells
        
        if ~cluInd(j);continue;end
        
        currSTimesBin = round(tempST{j} * ephysObj.params('APFs')) + sampOffset; % back to samples using same Fs as loadSpikesNPix; round, as ceil can be 1 off from float error
        currChannels  = max([1,ephysObj.cell_ID(j,3)-options.getnch]):min([nChanSorted, ephysObj.cell_ID(j,3)+options.getnch]); % sorted channel space; take care not to go <0 or > nChan
        
        % remove channels from list
        ind = ismember(currChannels,options.remchans);
        if any(ind)
            currChannels = currChannels(~ind);
            lhsAdd = sum(find(ind) <= options.getnch);
            if lhsAdd > 0;  currChannels = [currChannels(1)-lhsAdd:currChannels(1)-1, currChannels]; end
            rhsAdd = sum(find(ind) > options.getnch+1);
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
        if ~isempty(options.nwave)
            nWFs2Extract = min(length(currSTimesBin),options.nwave);
        else
            nWFs2Extract = length(currSTimesBin);
        end
        
        % drift corr. data only contains sorted channels; raw .ap.bin (single or cat) needs mapping back to raw channels
        if ~driftFlag
            currChannels = scanpix.npixUtils.mapChans(chanMap.connected,currChannels);
        end
        %
        currWave     = nan(nWFs2Extract,options.nsamp+1,2*options.getnch+1);
        ind2extract  = ceil(linspace(1,length(currSTimesBin),nWFs2Extract));

        % for peeling: index range of all sorted spikes whose template overlaps each snippet
        if peelFlag && nWFs2Extract > 0
            tSnip   = currSTimesBin(ind2extract);
            nearLo  = firstGE(tmpl.st, tSnip + wRel(1)   - tmpl.relT(end));
            nearHi  = firstGE(tmpl.st, tSnip + wRel(end) - tmpl.relT(1) + 1) - 1;
            ownClu  = ephysObj.cell_ID(j,1) - 1; % cell_ID is 1 based, spike_clusters.npy 0 based
        end

        c = 1;
        for k = ind2extract

            % skip spikes whose (filter) window would run over the file edges - leave as NaN
            if applyFilter; padSamp = [0.25*options.fs, 0.25*options.fs]; else; padSamp = [nSampPrePeak, options.nsamp-nSampPrePeak]; end
            if currSTimesBin(k)-padSamp(1) < 1 || currSTimesBin(k)+padSamp(2) > nSamp
                c = c+1;
                continue
            end

            if applyFilter
                startIdx = max([currSTimesBin(k)-0.25*options.fs,1]);
                endIdx   = min([currSTimesBin(k)+0.25*options.fs,nSamp]);
                currData = double(mmf.Data.x(:,startIdx:endIdx)'); %* winv;
                currData = HPfilter(currData, currChannels, find(logical(chanMap.connected)), options.fs, options.fhigh, options.flow, options.forder);
                currData = currData(0.25*options.fs-nSampPrePeak+1:0.25*options.fs+options.nsamp-nSampPrePeak+1,currChannels);
            else
                startIdx = max([currSTimesBin(k)-nSampPrePeak,1]);
                endIdx   = min([currSTimesBin(k)+options.nsamp-nSampPrePeak,nSamp]);
                currData = double(mmf.Data.x(:,startIdx:endIdx)');
                if driftFlag
                    currData = currData * winv(:,currChannels); % unwhiten (only output channels needed)
                else
                    currData = currData(:,currChannels);
                end
            end

            % subtract amplitude scaled templates of all other clusters' spikes overlapping this snippet
            if peelFlag
                currData = currData - reconstructOthers(tmpl, nearLo(c):nearHi(c), currSTimesBin(k), wRel, currChannels, ownClu);
            end
            currWave(c,:,1:size(currData,2)) = currData .* options.bitres ./ options.gain .* 1e6; %% uV conversion
            %
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
    if options.save
        waveforms = [tmpWaveforms tmpChannels];
        save(fullfile(ephysObj.dataPath{i},'waveforms.mat'),'waveforms','-v7.3');
        %
        path2ChanMap = dir(fullfile(fileparts(ephysObj.dataPath{i}),'*kilosortChanMap.mat'));
        try
            chanMap  = load(fullfile(path2ChanMap.folder,path2ChanMap.name));
            chansOrg = scanpix.npixUtils.mapChans(chanMap.connected,chanMap.chanMap);
            save(fullfile(ephysObj.dataPath{i},'waveForms_channelsRaw.mat'),'chansOrg','-v7.3');
        catch
            warning('scaNpix::ephys::extract_waveforms:Can''t find channel map that was used for kilosorting the data. Can''t generate ''waveForms_channelsRaw.mat''.');
        end
        %
        cellIDs = ephysObj.cell_ID;
        save(fullfile(ephysObj.dataPath{i},'waveForms_cluIDs.mat'),'cellIDs','-v7.3');
    end
    
end
close(hWait);

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
function filtData = HPfilter(data, channels, refChans, fs, fHigh, fLow, order)
% follows the Kilosort filter (function 'gpufilter' from the kilosort repo): mean subtraction, CAR by median, zero
% phase butterworth. Cut-offs are parameters so they can be matched to the preprocessing of the sorted data
% (e.g. CatGT -apfilter=butter,12,300,9000 -gblcar). fLow = [] or >= fs/2 -> high pass only (kilosort default)
% refChans: channels used for CAR - connected (neural) channels only, i.e. excluding ref and sync channels

% set up the parameters of the filter
if ~isempty(fLow) && fLow < fs/2
    [b1, a1] = butter(order, [fHigh/fs,fLow/fs]*2, 'bandpass');
else
    [b1, a1] = butter(order, fHigh/fs*2, 'high');
end

% subtract the mean from each channel
data = data - mean(data, 1); % subtract mean of each channel

% CAR, common average referencing by median across neural channels
data = data - median(data(:,refChans), 2); % subtract median across channels

% next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)
data(:,channels) = filter(b1, a1, data(:,channels)); % causal forward filter
data(:,channels) = flipud(data(:,channels)); % reverse time
data(:,channels) = filter(b1, a1, data(:,channels)); % causal forward filter again
filtData = flipud(data); % reverse time back

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tmpl = loadTemplates(ksDir, mmf, winv, nSamp, tAlign)
% load sorting output of concatenated data and prepare unwhitened, scaled templates for peeling
% Approach as in neuropixel-utils (KilosortDataset.reconstructRawSnippetsFromTemplates), but with KS3 output
% (rezToPhy2) amplitude x template does not map onto the data with one global factor (KS fits data/scaleproc and
% templates are re-estimated after amplitudes are fitted), so each template's scale is calibrated on the data here.

nCalib = 50; % n spikes/template used for scale calibration

T          = double(readNPY(fullfile(ksDir,'templates.npy')));   % nTemplates x nTimePoints x nChanSorted (whitened space)
tmpl.amp   = double(readNPY(fullfile(ksDir,'amplitudes.npy')));
tmpl.st    = double(readNPY(fullfile(ksDir,'spike_times.npy'))); % samples in cat file
tmpl.id    = double(readNPY(fullfile(ksDir,'spike_templates.npy'))) + 1;
tmpl.clu   = double(readNPY(fullfile(ksDir,'spike_clusters.npy'))); % 0 based, as phy
[tmpl.st, sortInd] = sort(tmpl.st);
tmpl.amp   = tmpl.amp(sortInd); tmpl.id = tmpl.id(sortInd); tmpl.clu = tmpl.clu(sortInd);
[nTmp, nT, nCh] = size(T);

% template sample aligned to spike time: trough of templates on their peak channel (KS pads templates so this
% is the same for all; checked against data for KS3: sample 41 of 82)
[~, pkCh] = max(squeeze(max(abs(T),[],2)),[],2);
if isempty(tAlign)
    trough = arrayfun(@(t) find(T(t,:,pkCh(t)) == min(T(t,:,pkCh(t))),1), (1:nTmp)');
    tAlign = median(trough(any(T(:,:),2)));
end
tmpl.relT = (1:nT) - tAlign; % template samples rel. to spike time

% calibrate scale of amp x template for each template against the (whitened) data on its 7 peak channels
k = nan(nTmp,1);
for t = 1:nTmp
    idx = find(tmpl.id == t);
    if numel(idx) < 10; continue; end
    idx  = idx(round(linspace(1,numel(idx),min(nCalib,numel(idx)))));
    rows = max(1,pkCh(t)-3):min(nCh,pkCh(t)+3);
    P    = reshape(T(t,:,rows), nT, numel(rows));
    [num, den] = deal(0);
    for s = idx'
        cols = tmpl.st(s) + tmpl.relT;
        if cols(1) < 1 || cols(end) > nSamp; continue; end
        D   = double(mmf.Data.x(rows, cols))';
        D   = D - mean(D(1:find(any(P,2),1),:),1); % baseline from zero padded part of template
        p   = tmpl.amp(s) * P;
        num = num + sum(D .* p,'all');
        den = den + sum(p.^2,'all');
    end
    k(t) = num/den;
end
k(isnan(k)) = median(k,'omitnan');
tmpl.scale  = k;

% scaled + unwhitened templates (same units as unwhitened data)
tmpl.Tu = reshape(reshape(T .* k, nTmp*nT, nCh) * winv, nTmp, nT, nCh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recon = reconstructOthers(tmpl, nearInd, t, wRel, chans, ownClu)
% sum of amplitude scaled templates of spikes (other than from ownClu) that overlap snippet around spike time t

nT    = numel(tmpl.relT);
tAl   = -tmpl.relT(1) + 1;   % template index of spike time
recon = zeros(numel(wRel), numel(chans));
for s = nearInd
    if tmpl.clu(s) == ownClu; continue; end
    ti    = wRel - (tmpl.st(s) - t) + tAl;   % template sample for each snippet sample
    valid = ti >= 1 & ti <= nT;
    recon(valid,:) = recon(valid,:) + tmpl.amp(s) * reshape(tmpl.Tu(tmpl.id(s), ti(valid), chans), nnz(valid), numel(chans));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = firstGE(x, a)
% vectorised binary search: for each a, index of first element in sorted x with x >= a (numel(x)+1 if none)

lo = ones(size(a));
hi = (numel(x)+1) * ones(size(a));
act = lo < hi;
while any(act)
    mid       = floor((lo(act) + hi(act)) / 2);
    goRight   = x(mid) < a(act);
    l = lo(act); h = hi(act);
    l(goRight)  = mid(goRight) + 1;
    h(~goRight) = mid(~goRight);
    lo(act) = l; hi(act) = h;
    act = lo < hi;
end
ind = lo;
end
