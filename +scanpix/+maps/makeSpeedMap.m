function speedMap = makeSpeedMap(spikeTimes,speed,trialDur,varargin)
%UNTITLED2 Summary of this function goes here


%% parse input
prms.posRate        = 50;    % Hz (50)
prms.minBinProp     = 0.005; % valid speed bins need to contain > prctg of samples of population (0.5%)
prms.speedBinSz     = 2;     % cm/s (2 cm/s)
prms.maxSpeed       = 40;    % cm/s (40 cm/s)
prms.smKernelLength = 0.25;  % in seconds (250ms)
prms.normaliseFR    = true;  % logical flag
prms.confInt        = 95;    % confidence interval


%% parse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(spikeTimes)
    spikeTimes = {spikeTimes};
end

% anon. fcn for 95% conf interval calculation - assuming normal dist which
% is prob. not the best way
CI = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x(:)))-1); %%

%% make speed map

% bin running speeds
[instSpeed,~,ind] = histcounts(speed,0:prms.speedBinSz:prms.maxSpeed);
validBins = instSpeed > prms.minBinProp * trialDur*prms.posRate; % in Kropf et al. it's >5% of all bins

speedMap = cell(length(spikeTimes),3);

for i = 1:length(spikeTimes)
    % inst firing rate
    instFRate = histcounts(ceil(spikeTimes{i}*prms.posRate),0:trialDur*prms.posRate) .* prms.posRate;
    % smooth
    kernel    = ones(1,ceil(prms.smKernelLength * prms.posRate)) ./ (prms.smKernelLength * prms.posRate); % Kropf uses kernel of 250ms
    instFRate = imfilter(instFRate,kernel,'replicate');
    
    % mean rate / speed bin 
    speedMap{i,1}(:,1) = prms.speedBinSz/2:prms.speedBinSz:prms.maxSpeed-prms.speedBinSz/2;
    speedMap{i,1}(:,2) = accumarray(ind(ind~=0),instFRate(ind~=0)',[length(instSpeed) 1],@mean);
    % 95 CI for means (a bit dense)
    speedMap{i,1}(:,3:4) = cell2mat( cellfun(@(x) CI(x,prms.confInt), accumarray(ind(ind~=0),instFRate(ind~=0)',[length(instSpeed) 1],@(x) {x}),'uni',0) );
    speedMap{i,1}        = speedMap{i,1}(validBins,:);
    % add Kropf et al. normalisation
    if prms.normaliseFR
        b = regress(speedMap{i}(:,2),[speedMap{i,1}(:,1) ones(length(speedMap{i,1}(:,1)), 1)]);
        yInt = speedMap{i}(1,1) - b(1)*prms.speedBinSz/2;
        speedMap{i,2} = (speedMap{i}(:,2) - yInt) ./ (b(1)*prms.maxSpeed);
        speedMap{i,3} = [b(1) yInt];
    end
end


end

