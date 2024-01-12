function [syncTTLs, missedSyncs] = loadSyncData(obj, trialIterator, varargin)
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
mode         = 'catgt';
lfpFs        = 2500;
syncChannel  = 385;
syncSigBit   = 10;
%
p = inputParser;
addParameter(p,'mode',     mode,        @ischar);
addParameter(p,'fs',       lfpFs,       @isscalar);
addParameter(p,'syncchan', syncChannel, @isscalar);
addParameter(p,'syncbit',  syncSigBit,  @isscalar);
%
parse(p,varargin{:});

%% process sync channel data
switch p.Results.mode
    case 'catgt'
        % first try looking for .txt file output from CatGT...
        syncTTLFile = dir(fullfile(obj.dataPath{trialIterator},'*TTL.txt'));

        if isempty(syncTTLFile)
            [fName,fPath] = uigetfile(fullfile(obj.dataPath{trialIterator},'*.txt'),'Select txt file with sync TTL times for current data set');
            cd(fPath);
        else
            fName = syncTTLFile.name;
            cd(syncTTLFile.folder);
        end

        fid = fopen(fName);
        syncTTLs = cell2mat(textscan(fid,'%f'));
        % syncTTLs = syncTTLs{1};
        fclose(fid);
    case 'lfp'
        cd(obj.dataPath{trialIterator});
        syncTTLs = loadSyncFromLFP(nChannels,p.Results.syncchan,p.Results.syncbit,p.Results.fs);
    otherwise
        error('scaNpix::ephys::loadSyncData:%s is not a valid method to load the sync data.', p.Results.mode)
end

% rarely there are some missing sync pulses in the npix stream
missedSyncs = [];
if ~obj.isConcat
    missedSyncs = find(diff(syncTTLs) > 1.5*1/obj.params('posFs'));
    if ~isempty(missedSyncs)
        missedPulses = [];
        missedSyncs(:,2) = floor((syncTTLs(missedSyncs+1) - syncTTLs(missedSyncs)) * obj.params('posFs'));
        missedSyncs(:,3) = syncTTLs(missedSyncs(:,1)) - syncTTLs(1);
        cs_NMissed = cumsum([0;missedSyncs(:,2)]);
        x = 1:length(syncTTLs);
        for i = 1:size(missedSyncs,1)
            missedPulses    = [missedPulses,(missedSyncs(i,1)+1:missedSyncs(i,1)+missedSyncs(i,2)) + cs_NMissed(i)];
            x(missedSyncs(i)+1:end) = x(missedSyncs(i,1)+1:end) + missedSyncs(i,2); % bump sample points for interpolation
        end
        interp_pulseT          = interp1(x, syncTTLs', missedPulses);
        temp                   = zeros(length(syncTTLs)+length(missedPulses),1);
        temp(missedPulses,1)   = interp_pulseT;
        temp(temp(:,1) == 0,1) = syncTTLs;
        syncTTLs               = temp;
        warning('scaNpix::ephys::loadSyncData: %i missing sync pulses across %i chunks in neuropixel datastream. Interpolated missing samples, but better go and check that',length(missedPulses),i);
    end
end

end

% wrapper to load sync data from lfp
function syncTTLs = loadSyncFromLFP(nChannels,syncChannel,syncSigBit,lfpFs)

try
    tmp = load('syncTTLs.mat');
    f = fieldnames(tmp);
    syncTTLs = tmp.(f{1}){2}; % only us
catch
    syncChData = scanpix.npixUtils.extractSyncChannel(cd, nChannels, syncChannel); % from cortex lab repository - loads sync channel from LFP file
    eventTimes = scanpix.npixUtils.spikeGLXdigitalParse(syncChData, lfpFs ); % from cortex lab repository - demultiplexes the 16 bit sync channel - SLOW!
    syncTTLs   = eventTimes{syncSigBit}; % {1} = all; {2} = ON; {3} = OFF ; times are in s
    save('syncTTLs.mat','syncTTLs'); % save to disk as a huge time sink
    syncTTLs = syncTTLs{2};
end

end

