function loadLFPNPix(obj, trialIterator)
% loadLFPNPix - load LFP data from neuropixel 1.0 files (.lf.bin)
% Data are kept as raw int16 to save memory/disk space. Multiply by
% obj.trialMetaData(trialIterator).lfpUVPerBit to convert to uV.
% Which channels get loaded is controlled by obj.lfpParams (see
% scanpix.helpers.defaultParamsLFP)
%
% package: scanpix.npixUtils
%
% Syntax:  loadLFPNPix(obj, trialIterator)
%
% Inputs:
%    obj           - ephys class object ('npix')
%    trialIterator - numeric index for trial to be loaded
%
% Outputs:
%    obj.lfpData.lfp{trialIterator}               - int16 [nChannels x nSamples]
%    obj.lfpData.lfpChans{trialIterator}          - [probe channel (1-based), depth from tip in um]; sorted by depth (ascending, like cell_ID)
%    obj.trialMetaData(trialIterator).lfpFs       - sample rate of stored LFP
%    obj.trialMetaData(trialIterator).lfpUVPerBit - [1 x nChannels] conversion factor to uV
%
%    1st LFP sample is the one closest to trial start (i.e. assume sample 1 = t0, as for DACQ)
%
% See also: scanpix.helpers.defaultParamsLFP
%
% LM 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialIterator (1,1) {mustBeNumeric}
end

%% params
prms = scanpix.helpers.defaultParamsLFP;
f    = fieldnames(obj.lfpParams);
for i = 1:length(f)
    prms.(f{i}) = obj.lfpParams.(f{i}); % user settings take precedence
end
if ~(isscalar(prms.downsampleFactor) && prms.downsampleFactor >= 1 && mod(prms.downsampleFactor,1) == 0)
    error('scaNpix::ephys::loadLFPNPix:''downsampleFactor'' needs to be a positive integer.');
end

fprintf('Loading LFP Data for %s .......... ', obj.trialNames{trialIterator});

%% find files
lfFile = fullfile(obj.dataPath{trialIterator}, [char(obj.trialNames{trialIterator}) '.lf.bin']);
if ~isfile(lfFile)
    lfFileInfo = dir(fullfile(obj.dataPath{trialIterator},'*.lf.bin'));
    if isscalar(lfFileInfo)
        lfFile = fullfile(lfFileInfo.folder, lfFileInfo.name);
    else
        [fName, fPath] = uigetfile(fullfile(obj.dataPath{trialIterator},'*.lf.bin'),['Select .lf.bin file for ' char(obj.trialNames{trialIterator})]);
        if isnumeric(fName)
            warning('scaNpix::ephys::loadLFPNPix:No LFP file selected. LFP for %s not loaded.', obj.trialNames{trialIterator});
            return
        end
        lfFile = fullfile(fPath, fName);
    end
end
metaFile = [lfFile(1:end-4) '.meta'];
if ~isfile(metaFile)
    error('scaNpix::ephys::loadLFPNPix:Can''t find meta file %s.', metaFile);
end

%% parse meta data
meta = readMeta(metaFile);

%%
if isfield(meta,'imDatPrb_type') && ~ismember(str2double(meta.imDatPrb_type),[0 1100 1300])
    error('scaNpix::ephys::loadLFPNPix:Probe type %s is not a neuropixel 1.0 probe. Only NP 1.0 has a separate LF stream.', meta.imDatPrb_type);
end
%
nChanFile  = str2double(meta.nSavedChans);
acqCounts  = str2num(meta.acqApLfSy); %#ok<ST2NM> % [AP LF SY] acquired
fileCounts = str2num(meta.snsApLfSy); %#ok<ST2NM> % [AP LF SY] saved in this file (AP = 0 for .lf.bin)
if fileCounts(2) == 0
    error('scaNpix::ephys::loadLFPNPix:No LF channels saved in %s.', lfFile);
end

%%

% file columns are ordered AP, LF, SY - LF columns follow any AP columns
lfCols = fileCounts(1) + (1:fileCounts(2));
% original (0-based) channel IDs of the columns in the binary file
if strcmp(meta.snsSaveChanSubset,'all')
    origChans = 0:nChanFile-1;
else
    origChans = str2num(meta.snsSaveChanSubset); %#ok<ST2NM> % e.g. '384:767,768'
end
lfIDs = origChans(lfCols);
if any(lfIDs >= acqCounts(1))
    lfIDs = lfIDs - acqCounts(1); % subset given in SpikeGLX's global numbering (LF = acqAP..acqAP+acqLF-1)
end
lfProbeChans = lfIDs + 1; % 1-based probe channel

%%
% gain per channel from imro table: (chan bank refID apGain lfGain [apHiPassFlt])
imroEntries = regexp(meta.imroTbl,'\(([^()]*)\)','tokens');
imroEntries = cellfun(@(x) str2num(x{1}), imroEntries(2:end),'uni',0); %#ok<ST2NM> % 1st entry is header
imro        = vertcat(imroEntries{:});
lfGain      = nan(acqCounts(2),1);
lfGain(imro(:,1)+1) = imro(:,5);
imAiRangeMax = 0.6;
if isfield(meta,'imAiRangeMax'); imAiRangeMax = str2double(meta.imAiRangeMax); end
imMaxInt     = 512; % 10 bit ADC on NP 1.0
if isfield(meta,'imMaxInt'); imMaxInt = str2double(meta.imMaxInt); end
uVPerBit     = imAiRangeMax ./ imMaxInt ./ lfGain .* 1e6;

% use nominal Fs - spike times are also converted with the nominal AP Fs
% and LF samples are locked to AP samples (1 LF sample = 12 AP samples),
% so this keeps spikes and LFP aligned
if isKey(obj.params,'lfpFs')
    fs = obj.params('lfpFs');
else
    fs = 2500;
end

%% channel selection
ycoords   = double(obj.chanMap(trialIterator).ycoords(:));
connected = logical(obj.chanMap(trialIterator).connected(:));
if numel(ycoords) ~= acqCounts(2)
    error('scaNpix::ephys::loadLFPNPix:Channel map has %d channels, but probe has %d. Channel map needs to cover the whole probe.', numel(ycoords), acqCounts(2));
end
availChans = intersect(find(connected), lfProbeChans); % excludes reference channels

if ~isempty(prms.chans)
    selChans = unique(prms.chans(:));
    if ~all(ismember(selChans, lfProbeChans))
        error('scaNpix::ephys::loadLFPNPix:Some of the requested channels in ''lfpParams.chans'' are not saved in %s.', lfFile);
    end
else
    if isnumeric(prms.depthRange)
        depthRange = sort(prms.depthRange(:))';
    elseif strcmpi(prms.depthRange,'cells') && ~isempty(obj.cell_ID)
        depthRange = [min(obj.cell_ID(:,2)) max(obj.cell_ID(:,2))];
    else
        if strcmpi(prms.depthRange,'cells')
            warning('scaNpix::ephys::loadLFPNPix:No spike data loaded, so can''t determine depth range from cells. Will span the whole probe instead.');
        end
        depthRange = [min(ycoords(availChans)) max(ycoords(availChans))];
    end
    % start at most dorsal position (depths are distance from probe tip)
    targetDepths = depthRange(2):-prms.chanSpacing:depthRange(1);
    selChans     = zeros(length(targetDepths),1);
    for i = 1:length(targetDepths)
        [~, ind]    = min(abs(ycoords(availChans) - targetDepths(i)));
        selChans(i) = availChans(ind);
    end
    selChans = unique(selChans);
end
% sort by depth (ascending, same as cell_ID)
[~, sortInd] = sortrows([ycoords(selChans) selChans]);
selChans     = selChans(sortInd);
[~, colInd]  = ismember(selChans, lfProbeChans);
readCols     = lfCols(colInd);

%% samples to read - trim to trial as done for spikes (i.e. 1st sync TTL to trial duration)
if isfield(obj.trialMetaData(trialIterator),'offSet') && ~isempty(obj.trialMetaData(trialIterator).offSet)
    offSet = obj.trialMetaData(trialIterator).offSet;
else
    offSet = 0;
    warning('scaNpix::ephys::loadLFPNPix:No sync offset available. LFP will not be aligned to trial start.');
end
fileInfo  = dir(lfFile);
nSampFile = fileInfo.bytes / 2 / nChanFile;
s0        = max(round(offSet * fs), 0);                                                         % 0-based index of 1st sample (closest to trial start)
s1        = min(floor((offSet + obj.trialMetaData(trialIterator).duration) * fs), nSampFile-1); % 0-based index of last sample
nSamp     = s1 - s0 + 1;

%% read data in chunks (file is interleaved: [nChanFile x nSamples])
fid     = fopen(lfFile,'r');
cleanUp = onCleanup(@() fclose(fid));
fseek(fid, s0 * nChanFile * 2, 'bof');

lfp       = zeros(length(selChans), nSamp, 'int16');
chunkSize = 1e5;
sampInd   = 1;
while sampInd <= nSamp
    nRead = min(chunkSize, nSamp - sampInd + 1);
    dat   = fread(fid, [nChanFile, nRead], '*int16');
    if size(dat,2) < nRead
        error('scaNpix::ephys::loadLFPNPix:Unexpected end of file in %s.', lfFile);
    end
    lfp(:,sampInd:sampInd+nRead-1) = dat(readCols,:);
    sampInd = sampInd + nRead;
end
clear cleanUp; % close file

%% downsample (optional) - same anti-aliasing filter as 'decimate', but keep 1st sample so timing stays simple
r = prms.downsampleFactor;
if r > 1
    [b, a]  = cheby1(8, 0.05, 0.8/r);
    lfpDS   = zeros(length(selChans), ceil(nSamp/r), 'int16');
    for i = 1:length(selChans)
        tmp        = filtfilt(b, a, double(lfp(i,:)));
        lfpDS(i,:) = int16(tmp(1:r:end));
    end
    lfp = lfpDS;
end

%% output
obj.lfpData(1).lfp{trialIterator}           = lfp;
obj.lfpData(1).lfpChans{trialIterator}      = [selChans ycoords(selChans)];
obj.trialMetaData(trialIterator).lfpFs       = fs / r;
obj.trialMetaData(trialIterator).lfpUVPerBit = uVPerBit(selChans)';

fprintf('  DONE! (%d channels)\n', length(selChans));
end

%%
function meta = readMeta(metaFile)
% parse spikeGLX meta file into struct (see also SGLXMetaToCoords_v2)
fid = fopen(metaFile, 'r');
C   = textscan(fid, '%[^=] = %[^\r\n]');
fclose(fid);
meta = struct();
for i = 1:length(C{1})
    tag = C{1}{i};
    if tag(1) == '~'
        tag = tag(2:end);
    end
    meta.(tag) = C{2}{i};
end
end
