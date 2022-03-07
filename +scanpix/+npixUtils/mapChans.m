function chansMapped = mapChans(connected,chans)
% mapChans - map channels from sorting to raw
% package: scanpix.npixUtils
%
% When you spike sort data with kilosort you will exclude at least the ref
% channel, but potentially more channels that might have been broken/noisy.
% Here we simply map channels back to the raw data channel ID
 %
%  Usage:   scanpix.npixUtils.mapChans( connected, chans ) 
%
%  Inputs:  
%           connected - logical array of connected channels for kilosort
%                       chan map (needs to be the full channel map)
%           chans     - numeric array of channels to be mapped
%
% LM 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chansMapped = chans;
remChan = find(~connected);
for i = 1:length(remChan)
    chansMapped(chansMapped>remChan(i)) = chansMapped(chansMapped>remChan(i)) + 1; 
end
end

