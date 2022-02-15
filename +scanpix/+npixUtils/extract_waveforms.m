function [npixObj,waveforms,channels] = extract_waveforms(npixObj,trialInd,varargin)
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

%%
p = parseParams(varargin); % see subfunction for parameter list

if nargin == 1 || isempty(trialInd)
    trialInd = 1;
end

%%
if strcmp(p.Results.mode,'drift')
    % load from drift corr file should be the default really whjen using
    % KS2.5 or 3
    if isempty(p.Results.path2drift)
        % prompt user to select file
        [fNameDrift,pathDrift] = uigetfile(fullfile(cd,'*.dat'),'Please Select Drift-corrected File');
        if isnumeric(pathDrift)
            warning('If you want to use the drift corrected data to extract waveforms, I need some info where that might be found on your disk. Too late now, but maybe you''ll do better later.');
            return
        end
        path2drift = [pathDrift fNameDrift];
    else
        path2drift = p.Results.path2drift;
    end
    % need a few details from concat log file
    fsLog   = dir(fullfile(fileparts(path2drift),'*logFile.tsv'));
    logFile = tdfread([fsLog.folder filesep fsLog.name],'tab');
    ind     = find(strcmp([npixObj.trialNames{trialInd} npixObj.fileType],cellstr(logFile.filename)));
    tempST  = cellfun(@(x) x + npixObj.trialMetaData(trialInd).offSet + sum(logFile.duration(2:ind-1)), npixObj.spikeData.spk_Times{trialInd},'uni',0); % add offset to spike times
%     tempST  = cellfun(@(x) x + npixObj.trialMetaData(trialInd).offSet, npixObj.spikeData.spk_Times{trialInd},'uni',0); % add offset to spike times

    [waveforms, channels] = getWaveforms(path2drift,tempST, npixObj.cell_ID(:,3), varargin);
elseif strcmp(p.Results.mode,'raw')
    % if you load from ap.bin raw - you should really HP filter and CAR this data before extracting waveforms
    tempST  = cellfun(@(x) x + npixObj.trialMetaData(trialInd).offSet, npixObj.spikeData.spk_Times{trialInd},'uni',0); % add offset to spike times
    if ~isempty(p.Results.clu)
        ind = ismember(npixObj.cell_ID(:,1),p.Results.clu);
    else
        ind = true(length(npixObj.cell_ID(:,1)),1);
    end
    [waveforms, channels] = getWaveforms(fullfile(npixObj.dataPath,[npixObj.trialNames{trialInd} npixObj.fileType]),tempST(ind), npixObj.cell_ID(ind,3), varargin);    
    
end
% assign to object directly as well
npixObj.spikeData.spk_waveforms{trialInd,1} = waveforms;
npixObj.spikeData.spk_waveforms{trialInd,2} = channels;
end

%% this subfunction grabs the actual waveform data
function [waveforms, channels] = getWaveforms(path2raw,spikeTimes, clu_Channel, varargin)
% getWaveforms - load waveform snippets from raw neuropixel data
% See other subfunction for parameter space
%
% Usage:
%       [waveforms, channels] = scanpix.npixUtils.getWaveforms(path2raw,spikeTimes, clu_Channel);
%       [waveforms, channels] = scanpix.npixUtils.getWaveforms(__, Name-Value comma separated list);
%
%
% Inputs:   path2raw    - path to raw data on disk
%           spikeTimes  - cell array of spike times for indivdual clusters  
%           clu_Channel - max amplitude channel ID for each cluster
%           varargin    - name-value: comma separated list of name-value pairs     
%
% Outputs:  waveforms   - array of waveforms for each cluster (format:)
%           channels    - array of channel IDs for each cluster
%
% LM 2020
%
%% TO DO:
% add HP filtering and CAR of raw file

%% PARAMS
p = parseParams(varargin{:});

if strcmp(p.Results.mode,'drift')
    unwhitenData = true; % need to unwhiten when using drift corr file!
   %%% Display some warning about bad param choice like using HP filter or CAR?
else
    unwhitenData = p.Results.unwhite; % this is prob redundant unless when you used KS2 and you want to load waveforms from ops.fproc
end

%% PREPROCESS
binFileStruct  = dir( path2raw );
if isempty(binFileStruct)
    error('Can''t find raw data mate...!')
end
dataTypeNBytes = numel(typecast(cast(0, p.Results.prec), 'uint8')); % determine number of bytes per sample
nSamp          = binFileStruct.bytes/(p.Results.nch*dataTypeNBytes);  % Number of samples per channel - quicker than reading from meta file

% we need to check in channel map if recording spanned multiple banks
chanMapFile    = dir( fullfile(binFileStruct.folder, '*ChanMap.mat') );
if isempty(chanMapFile)
    chanMapFName = scanpix.npixUtils.SGLXMetaToCoords_v2(binFileStruct.folder);
else
    chanMapFName = chanMapFile.name;
end
chanMap        = load(fullfile(binFileStruct.folder,chanMapFName));
bankBoundaries = find(abs(diff(chanMap.ycoords)) > p.Results.chanspace) + [0 1];

if unwhitenData
%     winv = readNPY(fullfile(binFileStruct.folder,'whitening_mat_inv.npy')); 
    try
        tmp  = load(fullfile(fileparts(path2raw),'whiteMat.mat')); % for KS3 we can't use the .npy anymore as it's just a diagonal matirx with the int16 scaling 
        winv = tmp.Wrot^-1;
    catch
        error('Can''t find whitening matrix! Either find that file or load from raw .ap.bin. Does that sound like a plan?');
    end
else
    winv = 1;  
end

nSampPrePeak = round(p.Results.prepeak*p.Results.nsamp);

%% load file
mmf = memmapfile(path2raw, 'Format', {p.Results.prec, [p.Results.nch nSamp], 'x'});
% remove common-mode noise
if p.Results.car
    %%%% THIS NEEDS TO MOVE INTO DEDICATED FUNCTION
    % remove ref (plus sync) and use only channels that were multi-
    % plexed together for median subtraction (every 24th in neuropix)
    %             allRefData                       = double(mmf.Data.x(:,currSTimesBin(j)-prms.nSamplesWf:currSTimesBin(j)+prms.nSamplesWf));
    %             allRefData(prms.chan2ignore,:)   = NaN;
    %             medVals                          = squeeze( nanmedian(reshape(allRefData(medInd,:)',2*prms.nSamplesWf+1,[],length(currChannels)),2) ); % medians per timepoint and channel
    %             % CAR
    %             currData                         = currData - medVals;
end

%% extract waveforms
hWait = waitbar(0);
[waveforms, channels] = deal(cell(length(spikeTimes),1));
for i = 1:length(spikeTimes) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    currSTimesBin = ceil(spikeTimes{i} * p.Results.fs);
    currChannels  = max([1,clu_Channel(i)-p.Results.getnch]):min([p.Results.nch, clu_Channel(i)+p.Results.getnch]); % take care not to go <0 or >prms.nCh
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
    
    currWave     = nan(nWFs2Extract,p.Results.nsamp+1,2*p.Results.getnch+1);
    ind2extract  = ceil(linspace(1,length(currSTimesBin),nWFs2Extract));
    c = 1;
    for j = ind2extract
        
        startIdx = max([currSTimesBin(j)-nSampPrePeak,1]);
        endIdx   = min([currSTimesBin(j)+p.Results.nsamp-nSampPrePeak,nSamp]);
        currData = double(mmf.Data.x(:,startIdx:endIdx)') * winv;
        currData = currData(:,currChannels);
        
        currWave(c,:,1:length(currChannels)) = currData .* p.Results.bitres ./ p.Results.gain .* 1e6; % uV conversion
        c = c+1;
    end
    waveforms{i} = currWave;
    channels{i}  = currChannels';
    
    waitbar( i/length(spikeTimes), hWait, ['cluster ' num2str(i) '/' num2str(length(spikeTimes))] );
end
close(hWait);

end

% subfunction to parse input
function p = parseParams(varargs)

defaultChanIgnore    = [];        % 192: reference channel; 385: sync channel - these should def be ignored
defaultPrecision     = 'int16';   % Data type of file
defaultNCh           = 385;       % Number of channels in file
defaultNChWaves      = 5;         % grab +/- this many channels around peak channel of cluster
defaultNWaves        = 250;       % this many waveforms/cluster (if [] we'll grab all)
defaultNSamplesWf    = 69;        % this many samples/AP (40 = 1.3ms)
defaultPropSampPre   = 0.4;       % relative amount of samples pre peak (0.375 @ 40 samples = 15 samples pre peak and 34 samples post peak
defaultFs            = 30000;     % sampleing rate
defaultGain          = 500;       % gain
defaultBitResolution = 1.2/2^10;  % in V/bit
defaultChanSpaceVert = 20;        % vertical channel spacing in um
defaultApplyCAR      = false;     % do common average referencing
defaultUnWhitedata   = false;     % unwhiten data - important if using KS2.5 or higher as drift corrected file is whitened  
mode                 = 'drift';
path2drift           = '';        % path to drift corr. file
cluIDs               = [];

p = inputParser;
addParameter(p,'remchans',   defaultChanIgnore);
addParameter(p,'prec',       defaultPrecision,@ischar);
addParameter(p,'nch',        defaultNCh,@isscalar);
addParameter(p,'getnch',     defaultNChWaves,@isscalar);
addParameter(p,'nwave',      defaultNWaves,@(x) isscalar(x) || isempty(x));
addParameter(p,'nsamp',      defaultNSamplesWf,@isscalar);
addParameter(p,'prepeak',    defaultPropSampPre);
addParameter(p,'fs',         defaultFs,@isscalar);
addParameter(p,'gain',       defaultGain,@isscalar);
addParameter(p,'bitres',     defaultBitResolution,@isscalar);
addParameter(p,'chanspace',  defaultChanSpaceVert,@isscalar);
addParameter(p,'car',        defaultApplyCAR);
addParameter(p,'unwhite',    defaultUnWhitedata);
addParameter(p,'mode',       mode,@ischar);
addParameter(p,'path2drift', path2drift,@ischar);
addParameter(p,'clu',        cluIDs);

parse(p,varargs{:});

end






