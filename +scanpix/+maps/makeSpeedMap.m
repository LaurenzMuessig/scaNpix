function speedMap = makeSpeedMap(spikeTimes,speed,trialDur,varargin)
%UNTITLED2 Summary of this function goes here


%% parse input
posRate = 50;
minBinProp = 0.005;
speedBinSz = 2; % cm/s
maxSpeed = 40;
smKernelLength = 0.25;

p = inputParser;
addOptional(p,'fs',posRate,@isscalar);
addParameter(p,'minprop',minBinProp,@isscalar);
addParameter(p,'binsz',speedBinSz,@isscalar);
addParameter(p,'maxspeed',maxSpeed,@isscalar);
addParameter(p,'smkernel',smKernelLength,@isscalar);
parse(p,varargin{:});

if ~iscell(spikeTimes)
    spikeTimes = {spikeTimes};
end

%% make speed map
[instSpeed,~,ind] = histcounts(speed,0:p.Results.binsz:p.Results.maxspeed);
validBins = instSpeed > p.Results.minprop * trialDur*p.Results.fs;

speedMap = cell(length(spikeTimes),1);
for i = 1:length(spikeTimes)
    instSpikeCount = histcounts(ceil(spikeTimes{i}*p.Results.fs),0:trialDur*p.Results.fs);
    % smooth
    kernel = ones(1,ceil(p.Results.smkernel * p.Results.fs)) ./ (p.Results.smkernel * p.Results.fs);
    instSpikeCount = conv(instSpikeCount,kernel,'same');
   
    % instSpikeCount2 = imfilter(instSpikeCount, kernel, 'replicate');
    
    speedMap{i} = accumarray(ind,instSpikeCount',[length(instSpeed) 1]) ./ (accumarray(ind,1,[length(instSpeed) 1]) / p.Results.fs);
    speedMap{i}(~validBins) = NaN;
    
    % add Kropf et al. normalisation!
    % [b,bint] = regress(y,X).

end



end

