% function [ampScaling, spkTemps] = loadSpikesNPix(obj, trialIterator, reloadFlag, loadAmps)
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
    obj.dataPathSort{trialIterator} = path2data_A;
else
    [path2data_A,path2data_B] = deal(obj.dataPathSort{trialIterator});
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
% spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - syncTTLs(1); % align to first TTL
spikeTimes         = double(spikeTimes) ./ obj.params('APFs') - obj.trialMetaData(trialIterator).offSet; % align to first TTL

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
% amps           = accumarray(clustIDs,amps,[max([unGoodClustIDs;good_clusts]) 1],@(x) {x});
% spkTemps       = accumarray(clustIDs,spkTemps,[max([unGoodClustIDs;good_clusts]) 1],@(x) {x});
% remove empty clusters - need to make sure not to remove cells that
% only fire in some trials of a trial sequence (we are assuming here that
% you clustered all data together and then split back into individual
% trials)
OtherClusters              = true(length(spikeTimesFin),1);
OtherClusters(good_clusts) = false;
indEmpty                   = cellfun('isempty',spikeTimesFin) & OtherClusters;
spikeTimesFin              = spikeTimesFin(~indEmpty);

%%%%%%%% I SHOULD CHECK WHETHER THIS WORKS FOR ALL CIRCUMSTANCES!!!!!!!
if ~loadFromPhy
    clu_Ch         = clu_Ch(~indEmpty);
    clu_Depth      = clu_Depth(~indEmpty);
end
% sort by depth
[clu_Depth, indSort] = sort(clu_Depth,'ascend'); % should be changed to descend to be sorted naturally
spikeTimesFin        = spikeTimesFin(indSort);

%% output
obj.spikeData(1).spk_Times{trialIterator} = spikeTimesFin;
if ~isempty(obj.posData.sampleT)
    obj.spikeData(1).sampleT{trialIterator}   = obj.spikeData(1).sampleT{trialIterator}(1:length(obj.posData(1).sampleT{trialIterator}));
else
    endIdxNPix                                = min( [ length(obj.spikeData(1).sampleT{trialIterator} ), find(obj.spikeData(1).sampleT{trialIterator} < obj.trialMetaData(trialIterator).duration,1,'last') + 1]);
    obj.spikeData(1).sampleT{trialIterator}   = obj.spikeData(1).sampleT{trialIterator}(1:endIdxNPix);
end

% what about reload?
if reloadFlag || trialIterator == 1 || isempty(obj.cell_ID)
    good_clusts    = good_clusts(indSort);
    cluLabel       = cluLabel(indSort);
    clu_Ch         = clu_Ch(indSort);
    % likely at least the ref channel will have been removed before sorting - this will map channel ID back to raw data
    clu_Ch_mapped  = scanpix.npixUtils.mapChans(obj.chanMap(trialIterator).connected,clu_Ch);
    
    obj.cell_ID    = [good_clusts, clu_Depth, clu_Ch clu_Ch_mapped];
    obj.cell_Label = cluLabel;
end

fprintf('  DONE!\n');
end

