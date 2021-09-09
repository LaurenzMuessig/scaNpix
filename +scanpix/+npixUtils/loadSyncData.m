function syncTTLs = loadSyncData(nFramesBonsai,BonsaiCorruptFlag,varargin)
% loadSyncData - load sync data to synchronise tracking data from Bonsai
% and neuropixel ephys data. We are assuming that the current directory
% contains the raw data and possibly output from CatGT
% package: scanpix.npixUtils
%
% Usage:
%       syncTTLs = scanpix.npixUtils.loadSyncData;
%       syncTTLs = scanpix.npixUtils.loadSyncData([],BonsaiCorruptFlag);
%       syncTTLs = scanpix.npixUtils.loadSyncData(nFramesBonsai,BonsaiCorruptFlag);
%       syncTTLs = scanpix.npixUtils.loadSyncData(__,prmsStruct);
%       syncTTLs = scanpix.npixUtils.loadSyncData(__, Name-Value comma separated list);
%
%
% Inputs:   nFramesBonsai     - n samples from Bonsai stream
%           BonsaiCorruptFlag - logical flag if Bonsai data is corrupt 
%           varargin          - prmsStruct: structure with parameter fields to be changed from defaults
%                             - name-value: comma separated list of name-value pairs     
%
% Outputs:  syncTTLs          - time in s for sync TTLs in neuropixel
%                               stream
%
% LM 2020


%% Params
prms.lfpFs         = 2500; % sampling rate for LFP, 2.5kHz is default
prms.nChannels     = 385;
prms.syncChannel   = 385;
prms.syncSigBit    = 10; % SMA input in IMEC board is bit #6 in SpikeGLX (counting down from 16)

% ---------------------------------------------------------------------------------- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %

if nargin < 2
    BonsaiCorruptFlag = false;
end

%% process sync channel data
% first try looking for .txt file output from CatGT...
syncTTLFile = dir('*TTL.txt');

if ~isempty(syncTTLFile)
    % best case - note that CatGT only registers complete TTLs, so
    % regularly drops last one
    fid = fopen(syncTTLFile.name);
    syncTTLs = textscan(fid,'%f');
    syncTTLs = syncTTLs{1};
    fclose(fid);    
else
    % if we didn't run CatGT (legacy data..)..
    % .. first try and load data
    syncTTLs = loadSyncFromLFP(prms);    
end

% in case we just wanted the sync data, we'll stop here and don't worry
% about Bonsai
if nargin == 0 || isempty(nFramesBonsai)
    return;
end

% this needs a bit more experimentation if we accounted for all eventualities
% in an acceptable manner, esp. for corrupt data
if ~BonsaiCorruptFlag && nFramesBonsai ~= length(syncTTLs)
    if nFramesBonsai - length(syncTTLs) == -1
        % this case isn't 100% clear as frame could be missing anywhere
        disp('Warning. Missmatch between n of pos samples and n of TTLs. -1 frame in pos data, so we assume the last frame in neuropix stream is incomplete and will be removed from neuropix data.');
        syncTTLs = syncTTLs(1:end-1);
    elseif nFramesBonsai - length(syncTTLs) == +1
        % this case should be clear and essentially no frame is missing!
        disp('Warning. Missmatch between n of pos samples and n of TTLs. +1 frame in pos data, so the last frame in neuropix stream is incomplete and will be removed from pos data.');
    elseif nFramesBonsai - length(syncTTLs) > 1
        % this case should happen when animal unplugs during recording
        fprintf('Warning. There are %i more frames in tracking stream compared to neuropix data - assuming that the animal unplugged in redcording. If not you are in trouble\n', nFramesBonsai-length(syncTTLs));
    else
        fprintf('Warning. Missmatch between n of pos samples and n of TTLs in neuropix data - %i vs. %i. Better go and check out why.\n', nFramesBonsai, length(syncTTLs));
    end
elseif BonsaiCorruptFlag
    % Note: It took me a while to figure out that FlyCap is sensitive to how
    % it's closed when switching between cameras or changing settings. This
    % results in the camera not sending it's metadata over (i.e. frame
    % count and time stamp). When we haven't got the framecount we can't
    % reconstruct which frames are missing. So by assuming we have an even
    % sampling interval we will introduce some jitter, so at some point the
    % 2 streams will be more and more out of sync (also depending on when 
    % the frames were skipped, the later the better). So for these cases we
    % should check the data, i.e. e.g. do we see an obvious drop in the 
    % spatial properties. 
    if nFramesBonsai == length(syncTTLs)
        % best case scenario - no frames are missing
        disp('Corrupt Point Grey MetaData Logging: No frame N mismatch with neuropix stream! All good. Phew...') 
    elseif nFramesBonsai - length(syncTTLs) == -1
        disp('Corrupt Point Grey MetaData Logging: -1 frame in position data. Assuming last frame is missing for Bonsai data...');
        syncTTLs = syncTTLs(1:end-1);
    elseif nFramesBonsai - length(syncTTLs) == +1
        disp('Corrupt Point Grey MetaData Logging: +1 frame in position data. Assuming last frame is missing in neuropix data...')
    else
        fprintf('Corrupt Point Grey MetaData Logging: Frame N mismatch is %i (Bonsai) vs. %i (Neuropix). We''ll assume even sampling between frames, but you should check carefully if n of mismatch is too big!\n', nFramesBonsai, length(syncTTLs));
    end
    
end

end

% wrapper to load sync data from lfp
function syncTTLs = loadSyncFromLFP(prms)

try
    tmp = load('syncTTLs.mat');
    f = fieldnames(tmp);
    syncTTLs = tmp.(f{1}){2}; % only us
catch
    syncChData = scanpix.npixUtils.extractSyncChannel(cd, prms.nChannels, prms.syncChannel); % from cortex lab repository - loads sync channel from LFP file
    eventTimes = scanpix.npixUtils.spikeGLXdigitalParse(syncChData, prms.lfpFs ); % from cortex lab repository - demultiplexes the 16 bit sync channel - SLOW!
    syncTTLs  = eventTimes{prms.syncSigBit}; % {1} = all; {2} = ON; {3} = OFF ; times are in s
    save('syncTTLs.mat','syncTTLs'); % save to disk as a huge time sink
    syncTTLs = syncTTLs{2};
end

end

