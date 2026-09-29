function prms = defaultParamsLFP
% defaultParamsLFP - Generate the default params structure for loading LFP
% data from neuropixel 1.0 recordings (.lf.bin stream)
% package: scanpix.helpers
%
% Syntax:  prms = scanpix.helpers.defaultParamsLFP
%
% Inputs:
%
% Outputs: prms - prms struct with default parameters
%
%          prms.chanSpacing      - spacing (in um) between loaded channels along the probe
%          prms.depthRange       - 'cells' - span from most dorsal to most ventral cell
%                                            (as per obj.cell_ID(:,2)); falls back to
%                                            'probe' if no spike data has been loaded yet
%                                  'probe' - span whole probe
%                                  [minDepth maxDepth] - depth range in um from probe tip
%          prms.chans            - numeric list of probe channels (1-based) to load;
%                                  overrides chanSpacing/depthRange if not empty
%          prms.downsampleFactor - integer factor to downsample LFP by (1 = no downsampling)
%
% See also: scanpix.npixUtils.loadLFPNPix
%
% LM 2026

prms.chanSpacing      = 100;     % in um
prms.depthRange       = 'cells'; % 'cells', 'probe' or [minDepth maxDepth] in um
prms.chans            = [];      % explicit channel list (1-based) - overrides the above
prms.downsampleFactor = 1;       % 1 = keep native Fs (2.5kHz)

end
