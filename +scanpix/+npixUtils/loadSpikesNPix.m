function loadSpikesNPix(obj, trialIterator, reloadFlag)
% loadSpikes - load spike data from neuropixel files
% We will just load spike times
%
% Syntax:  loadSpikes(obj, trialIterator)
%
% Inputs:
%    obj           - ephys class object ('npix')
%    trialIterator - numeric index for trial to be loaded
%
% Outputs:
%
% See also: 
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To-Do

% initialise
fprintf('Loading Neural Data for %s .......... ', obj.trialNames{trialIterator});

%% grab sync channel data

cd(obj.dataPath{trialIterator});
if reloadFlag
    syncTTLs = obj.trialMetaData(trialIterator).offSet;
elseif ~isfield(obj.trialMetaData,'BonsaiCorruptFlag')
    syncTTLs = scanpix.npixUtils.loadSyncData();
else
    syncTTLs = scanpix.npixUtils.loadSyncData(length(obj.posData.sampleT{trialIterator}),obj.trialMetaData(trialIterator).BonsaiCorruptFlag);
end

% in 1 single trial there are a bunch of sync pulses missing in middle of
% trial - not sure what happened there. hopefully it's a one off 
missedSyncs = find(diff(syncTTLs) > 1.5*1/obj.params('posFs'));
if ~isempty(missedSyncs)
    nMissedPulses = floor((syncTTLs(missedSyncs+1) - syncTTLs(missedSyncs)) * obj.params('posFs'));
    missedPulses  = missedSyncs+1:missedSyncs+nMissedPulses;
    interp_pulseT          = interp1([1:missedSyncs,missedSyncs+nMissedPulses+1:length(syncTTLs)+nMissedPulses], syncTTLs', missedPulses);
    temp                   = zeros(length(syncTTLs)+nMissedPulses,1);
    temp(missedPulses,1)   = interp_pulseT;
    temp(temp(:,1) == 0,1) = syncTTLs;
    syncTTLs               = temp;
    warning('scaNpix::ephys::loadSpikes:Missing sync pulses in neuropixel datastream. Interpolated missing samples, but better go and check that');
end


%             % decide what to load - phy or kilosort
%             if obj.params('loadFromPhy') && ~reloadFlag
%                 if exist(fullfile(obj.dataPath{trialIterator},'cluster_info.tsv'),'file') == 2
%                     loadFromPhy = true;
%                 else
%                     warning('scaNpix::npix::loadSpikes:Can''t find ''cluster_info'' from phy output. Will try using kilosort data instead!');
%                     loadFromPhy = false;
%                 end
%             else
%                 loadFromPhy = false;
%             end
%% load sorting resuts
% deal with directories in case of reload or KS output is in different directory to raw data
if reloadFlag || exist(fullfile(obj.dataPath{trialIterator},'spike_times.npy'),'file') ~= 2
    % for a relaod we always request user interaction
    [fName, path2data_A] = uigetfile(fullfile(obj.dataPath{trialIterator},'*.npy'),['Select spike_times.npy for ' obj.trialNames{trialIterator}]);
    if isnumeric(fName)
        warning('Cannot continue reloading spike sorting results. Your journey ends here young padawan!');
        return
    end
    path2data_B = path2data_A;
else
    [path2data_A,path2data_B] = deal(obj.dataPath{trialIterator});
end

% decide what to load - phy or kilosort
if obj.params('loadFromPhy')
    if exist(fullfile(path2data_A,'cluster_info.tsv'),'file') == 2
        loadFromPhy = true;
    else
        warning('scaNpix::ephys::loadSpikes:Can''t find ''cluster_info'' from phy output. Will try using kilosort data instead!');
        loadFromPhy = false;
    end
else
    loadFromPhy = false;
end

% try and load cluster IDs from back up folder if data was PHY'd as this will prevent issues further down in case you curated the data (and merged/split at least one cluster)
if ~loadFromPhy && exist(fullfile(obj.dataPath{trialIterator},'cluster_info.tsv'),'file') == 2
    
    if isfolder(fullfile(obj.dataPath{trialIterator},'backUpFiles'))
        path2data_A = fullfile(path2data_A,'backUpFiles');
    else
        warning('scaNpix::ephys::loadSpikes:Can''t find folder with BackUp files! If you merged/split clusters in PHY, loading raw kilosort results will fail shortly...');
    end
end

% load spike times
spikeTimes         = readNPY(fullfile(path2data_A,'spike_times.npy'));
spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
obj.trialMetaData(trialIterator).offSet = syncTTLs(1); % this is offset of first TTL from trial start - need  a record for anything related to raw data
% load cluster IDs
clustIDs           = readNPY(fullfile(path2data_A,'spike_clusters.npy')) + 1; % 0 based index

if loadFromPhy
    % phy curated output - only load clusters labelled good
    clu_info       = tdfread(fullfile(path2data_A,'cluster_info.tsv'),'tab');
    goodLabel      = strcmp(cellstr(clu_info.group),'good');
    try
        good_clusts    = clu_info.id(goodLabel) + 1;
    catch
        good_clusts    = clu_info.cluster_id(goodLabel) + 1;
    end
    clu_Depth      = clu_info.depth(goodLabel);  % I think this amd the following seem to not give correct results in somce cases?? Phy Issue??
    clu_Ch         = clu_info.ch(goodLabel) + 1; % this is 0 based
    cluLabel       = string(clu_info.group);
    cluLabel       = cluLabel(goodLabel);
else
    ks_labels = tdfread(fullfile(path2data_A,'cluster_KSLabel.tsv'),'tab'); % we need the raw kilosort output here before you ran Phy as this file gets overwritten! - maybe point to backup folder?
    % raw kilosort output - we'll take everything in this case and can
    % filter later if needed
    cluLabel       = string(ks_labels.KSLabel); % this is either 'good' or 'mua'
    %     ind_good       = strcmp(cluLabel,'good') | strcmp(cluLabel,'mua');
    good_clusts    = ks_labels.cluster_id + 1; % ID for clusters also 0 based
    % get depth estimate for each cluster based on template ampl
    % distribution across probe
    templates          = readNPY(fullfile(path2data_B,'templates.npy'));
    Winv               = readNPY(fullfile(path2data_B,'whitening_mat_inv.npy'));
    chanPos            = readNPY(fullfile(path2data_B,'channel_positions.npy'));
    chanMapKS          = double(readNPY(fullfile(path2data_B,'channel_map.npy'))) + 1;  %
    [clu_Depth,clu_Ch] = scanpix.npixUtils.getCluChDepthFromTemplates(templates, Winv, [chanMapKS(:) chanPos(:,2)]);
end

% now we need to remove bad clusters and spike times outside trial
unClustIDs     = unique(clustIDs);
% index for 'good' clusters
unGoodClustIDs = unClustIDs(ismember(unClustIDs,good_clusts)); % remove 'mua' or 'noise' clusters from list in case we deal with phy output
[~,indGood]    = ismember(clustIDs,unGoodClustIDs); % only keep these
% index for spike times within ]1st TTL, trialduration]
indSTimes      = spikeTimes > 0 & spikeTimes <= obj.trialMetaData(trialIterator).duration;
% final index of spike times we want to keep
ind2keep       = indGood ~= 0 & indSTimes;

% remove from data
clustIDs       = clustIDs(ind2keep);
spikeTimes     = spikeTimes(ind2keep);

% sort clusters, so accumarray output is sorted
[clustIDs, sortInd] = sort(clustIDs);
spikeTimes          = spikeTimes(sortInd);

% reformat into more convenient form
spikeTimesFin  = accumarray(clustIDs,spikeTimes,[max([unGoodClustIDs;good_clusts]) 1],@(x) {x});
% remove empty clusters - need to make sure not to remove cells that
% only fire in some trials of a trial sequence (we are assuming here that
% you clustered all data together and then split back into individual
% trials)
OtherClusters    = true(length(spikeTimesFin),1);
OtherClusters(good_clusts) = false;
indEmpty       = cellfun('isempty',spikeTimesFin) & OtherClusters;
spikeTimesFin  = spikeTimesFin(~indEmpty);

%%%%%%%% I SHOULD CHECK WHETHER THIS WORKS FOR ALL CIRCUMSTANCES!!!!!!!
if ~loadFromPhy
    clu_Ch         = clu_Ch(~indEmpty);
    clu_Depth      = clu_Depth(~indEmpty);
end
% sort by depth
[clu_Depth, indSort] = sort(clu_Depth,'ascend'); % should be changed to descent to be sorted naturally
spikeTimesFin        = spikeTimesFin(indSort);
%% output
obj.spikeData(1).spk_Times{trialIterator} = spikeTimesFin;
if ~reloadFlag
    endIdxNPix                                = min( [ length(syncTTLs), find(syncTTLs < obj.trialMetaData(trialIterator).duration + syncTTLs(1),1,'last') + 1]);
    obj.spikeData(1).sampleT{trialIterator}   = syncTTLs(1:endIdxNPix) - syncTTLs(1);
else
    endIdxNPix = size(obj.posData.XYraw{trialIterator},1);
end
% we'll only need to grab this once - this assumes data was clustered together, but wouldn't make much sense otherwise?
% what about reload?
if reloadFlag || trialIterator == 1
    good_clusts          = good_clusts(indSort);
    cluLabel             = cluLabel(indSort);
    clu_Ch               = clu_Ch(indSort);
    % likely at least the ref channel will have been removed before sorting - this will map channel ID back to raw data
    clu_Ch_mapped = scanpix.npixUtils.mapChans(obj.chanMap(trialIterator).connected,clu_Ch);
    
    obj.cell_ID    = [good_clusts, clu_Depth, clu_Ch clu_Ch_mapped];
    obj.cell_Label = cluLabel;
end

% also keep only tracking data within trial times
if ~isempty(obj.posData.XYraw) && ~isempty(obj.posData.XYraw{trialIterator})
    obj.posData.XYraw{trialIterator}     = obj.posData.XYraw{trialIterator}(1:endIdxNPix,:);
    obj.posData.XY{trialIterator}        = obj.posData.XY{trialIterator}(1:endIdxNPix,:);
    obj.posData.speed{trialIterator}     = obj.posData.speed{trialIterator}(1:endIdxNPix,:);
    obj.posData.direction{trialIterator} = obj.posData.direction{trialIterator}(1:endIdxNPix,:);
    obj.posData.sampleT{trialIterator}   = obj.posData.sampleT{trialIterator}(1:endIdxNPix,:);
end

fprintf('  DONE!\n');
end

