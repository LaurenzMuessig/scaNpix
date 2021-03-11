function preallocEmpty(obj,noLoadProps,loadProps)
% preallocEmpty - preallocate empty properties in dacq or npix class object 
% (e.g. from things that weren't loaded) that have a trial based structure
% (i.e. all raw data like positions, spiketimes etc.)
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.preallocEmpty(obj,noLoadProps)
%       scanpix.helpers.preallocEmpty(obj,noLoadProps,loadProps)
%
% Inputs:
%    noLoadProps - true/false - flag to do properties that won't be loaded
%                  through normal routines
%    loadProps   - cell array; any combination of {'pos','spikes','lfp'}
%
% Outputs:
%
% see also:
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    loadProps   = {};
end

% these will never get loaded automatically by normal loading routine
if noLoadProps
    obj.posData.linXY = cell(1,length(obj.trialNames));
    obj.maps          = structfun(@(x) cell(1,length(obj.trialNames)),obj.maps,'uni',0);
    obj.linMaps       = structfun(@(x) cell(1,length(obj.trialNames)),obj.linMaps,'uni',0);
end

% this could actually have a hard coded list instead of supplying as input?
if ~isempty(loadProps)
    for i = 1:length(loadProps)
        f = fieldnames(obj.(loadProps{i}));
        for j = 1:length(f)
            if isempty(obj.(loadProps{i}).(f{j}))
                obj.(loadProps{i}).(f{j}) = cell(1,length(obj.trialNames));
            end
        end
    end
end
end
