function [templateDepths, maxChan] = getCluChDepthFromTemplates(templates, wInv, chanPos)
% determine depth of template as well as closest channel. Adopted from
% 'templatePositionsAmplitudes' by N. Steinmetz (https://github.com/cortex-lab/spikes)
% package: scanpix.npixUtils
%
% Usage:
%       [templateDepths, maxChan] = scanpix.npixUtils.getCluChDepthFromTemplates(templates, wInv, chanPos);
%
% Inputs:   templates       - templates array from kilosort output ('templates.npy')
%           wInv            - inverse of whitening matrix ('whitening_mat_inv.npy') 
%           chanPos         - array with channel ID and channel position (['channel_map.npy' 'channel_positions.npy'(:,2)])
%
% Outputs:  templateDepths  - depth template for given clusters
%           maxChan         - closest channel ID for COM of given clusters
%
% LM 2020

%% TO DO: FIX ISSUE WITH KS2_5 ChanMap - First entry depth=NaN?

% unwhiten all the templates
tempsUnW = zeros(size(templates));
for t = 1:size(templates,1)
    tempsUnW(t,:,:) = squeeze(templates(t,:,:)) * wInv;
end

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel
tempAmpsUnscaled = max(tempChanAmps,[],2);

% need to zero-out the potentially-many low values on distant channels ...
threshVals = tempAmpsUnscaled*0.3; 
tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;

% ... in order to compute the depth as a center of mass
templateDepths = sum(bsxfun(@times,tempChanAmps,chanPos(:,2)'),2)./sum(tempChanAmps,2);

% closest channel to template depth 
[~, maxChanInd] = min(abs(templateDepths-chanPos(:,2)'),[],2);
maxChan = chanPos(maxChanInd,1);

end

