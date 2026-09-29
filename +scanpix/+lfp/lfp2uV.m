function [lfpUV, t, chans] = lfp2uV(obj, trialInd, options)
% lfp2uV - convert neuropixel LFP stored as raw int16 in object into uV
% package: scanpix.lfp
%
% Usage:
%       [lfpUV, t, chans] = scanpix.lfp.lfp2uV(obj, trialInd)
%       [lfpUV, t, chans] = scanpix.lfp.lfp2uV(obj, trialInd, 'chans', [41 61])
%       [lfpUV, t, chans] = scanpix.lfp.lfp2uV(obj, trialInd, 'precision', 'single')
%
% Inputs:   obj       - ephys class object ('npix') with LFP loaded
%           trialInd  - numeric index of trial
%
% Options:
%   'chans'      []        - probe channels (1-based, as in obj.lfpData.lfpChans{trialInd}(:,1))
%                            to convert; empty = all loaded channels
%   'precision'  'double'  - 'double' or 'single'
%
% Outputs:  lfpUV - LFP in uV [nChannels x nSamples]
%           t     - time (s) of each sample relative to trial start [1 x nSamples]
%           chans - [probe channel, depth in um] for each row of lfpUV
%
% See also: scanpix.npixUtils.loadLFPNPix
%
% LM 2026
%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeInteger, mustBePositive}
    options.chans {mustBeNumeric} = []
    options.precision (1,:) char {mustBeMember(options.precision,{'double','single'})} = 'double'
end

if ~strcmp(obj.type,'npix')
    error('scaNpix::lfp::lfp2uV:Only works for npix objects. DACQ LFP is already stored in volts.');
end
if trialInd > length(obj.lfpData.lfp) || isempty(obj.lfpData.lfp{trialInd})
    error('scaNpix::lfp::lfp2uV:No LFP loaded for trial %d.', trialInd);
end

chans = obj.lfpData.lfpChans{trialInd};
if isempty(options.chans)
    rowInd = 1:size(chans,1);
else
    [isLoaded, rowInd] = ismember(options.chans(:), chans(:,1));
    if ~all(isLoaded)
        error('scaNpix::lfp::lfp2uV:Channel(s) %s not loaded for trial %d.', num2str(options.chans(~isLoaded)'), trialInd);
    end
end

scale = cast(obj.trialMetaData(trialInd).lfpUVPerBit(rowInd), options.precision);
lfpUV = cast(obj.lfpData.lfp{trialInd}(rowInd,:), options.precision) .* scale(:);
chans = chans(rowInd,:);

t = (0:size(lfpUV,2)-1) ./ obj.trialMetaData(trialInd).lfpFs;

end
