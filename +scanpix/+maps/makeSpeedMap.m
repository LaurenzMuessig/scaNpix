function speedMaps = makeSpeedMap(obj, trialInd, cellInd, options)
% makeSpeedMap - generate running speed firing rate maps from spike times 
% and running speed of animal as in Kropf et al. (https://www.nature.com/articles/nature14622)
% Note that output is a bit dense - check what data is where in cell array
%
% package: scanpix.maps
%
% Syntax:
%       scanpix.maps.makeSpeedMap(spikeTimes,speed,trialDur)
%       scanpix.maps.makeSpeedMap(__, prmsStruct )
%       scanpix.maps.makeSpeedMap(__, Name-Value comma separated list )
%
% Inputs:
%    spkTimes    - 1xnCell cell array of spike times in seconds
%    speed       - (optional) nPosSamplesx1 array of running speeds in cm/s 
%    trialDur    - trial duration in seconds             
%    varargin    - prmsStruct: structure with parameter fields to be changed from defaults
%                - name-value: comma separated list of name-value pairs
%
% Outputs:
%   speedMaps    - speed map; cell array(nCells x 3) 
%                - dim1: (speedBins,speedMap,+ConfInterval,-ConfInterval)
%                - dim2: (normalised speedMap)
%                - dim3: (slope,Yintercept) of linear regression of speedMap
%
% see also: scanpix.maps.addMaps;
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
    cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);
    options.minSpeed (1,1) {mustBeNumeric} = 2;
    options.smooth (1,1) {mustBeNumericOrLogical} = false;
end


%%
spkTimes = obj.spikeData.spk_Times{trialInd}(cellInd);
if islogical(cellInd)
    outInd = find(cellInd);
else
    outInd = cellInd;
end

%%
% anon. fcn for 95% conf interval calculation - assuming normal dist which
% is prob. not the best way
CI = @(x,p) std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x(:)))-1);

%% make speed map

% bin running speeds
[speedHist,bins,ind]   = histcounts(obj.posData.speed{trialInd},linspace(0,obj.mapParams.speed.maxSpeed,ceil(obj.mapParams.speed.maxSpeed/obj.mapParams.speed.binSizeSpeed)+1));
lastValidBin           = find(speedHist > obj.mapParams.speed.minBinProp * sum(speedHist),1,'last'); % find last valid bin for data (in Kropf et al. it's >0.5% of all bins)
%
maxSpeed               = bins(lastValidBin);
ind(ind>=lastValidBin) = 0;
%
occCounts              = accumarray(ind(ind~=0), 1, [lastValidBin-1 1])./ obj.trialMetaData(trialInd).posFs;
%
speedMaps              = cell(length(obj.spikeData.spk_Times{trialInd}),2);
speedBins              = obj.mapParams.speed.binSizeSpeed/2:obj.mapParams.speed.binSizeSpeed:maxSpeed-obj.mapParams.speed.binSizeSpeed/2;
% smoothing kernel inst firing rate
kernel                 = fspecial('gaussian',[ceil(0.25 * obj.trialMetaData(trialInd).posFs) 1],2); % Kropf et al. use kernel of 250ms

%
if obj.mapParams.speed.showWaitBar; hWait = waitbar(0); end

for i = 1:length(spkTimes)

    %
    if isempty(spkTimes{i})
        speedMaps{outInd(i),1}      = nan(length(speedBins),3);
        speedMaps{outInd(i),1}(:,1) = speedBins;
        speedMaps{outInd(i),2}      = NaN;  
        continue; 
    end

    % spike counts per speed bin
    instSpikeCount                  = accumarray(ceil(spkTimes{i} .* obj.trialMetaData(trialInd).posFs),1,size(obj.posData.speed{trialInd}));
    spikeCounts                     = accumarray(ind(ind~=0),instSpikeCount(ind~=0)',[lastValidBin-1 1]);

    % speed tuning curve
    speedMaps{outInd(i),1}(:,1)     = speedBins;
    speedMaps{outInd(i),1}(:,2)     = spikeCounts ./ occCounts;
    if options.smooth
        speedMaps{outInd(i),1}(:,2) = smoothdata(speedMaps{outInd(i),1}(:,2),1, "sgolay", 13, "Degree", 3);   
    end
    % confidence interval
    confInt                         = cell2mat( cellfun(@(x) CI(x,obj.mapParams.speed.confInt), accumarray(ind(ind~=0),instSpikeCount(ind~=0)',[lastValidBin-1 1],@(x) {x}),'uni',0) ) .* obj.trialMetaData(trialInd).posFs;
    speedMaps{outInd(i),1}(:,3)     = confInt(:,2);
    % speedMaps{outInd(i),1}(:,4)     = occCounts;
   
    % speed score (r speed vs instantaneous firing rate)
    instFRate                       = imfilter(instSpikeCount,kernel,'replicate');
    valSpeedInd                     = ~isnan(obj.posData.speed{trialInd}) & obj.posData.speed{trialInd} > options.minSpeed; % exclude periods of immobility
    speedMaps{outInd(i),2}          = corr(obj.posData.speed{trialInd}(valSpeedInd),instFRate(valSpeedInd));   
    %
    if obj.mapParams.speed.showWaitBar; waitbar(i/length(spkTimes),hWait,sprintf('Making those Speed Maps... %i/%i done.',i,length(spkTimes))); end
end

if obj.mapParams.speed.showWaitBar; close(hWait); end

end

 % linear regression
    % speedMaps{i,1}(:,4)     = NaN;
    % speedMaps{i,1}(1:2,4)   = regress(speedMaps{i,1}(:,2),[ones(length(speedMaps{i,1}(:,1)), 1) speedMaps{i,1}(:,1)]);

    % add Kropf et al. normalisation - probably pointless as can just compute from regression results
    % if obj.mapParams.speed.normaliseFR
    %     % speedMaps{i,1}(1:2,4) = (speedMaps{i,1}(:,2) - b(1) ) ./ (b(2) * obj.mapParams.speed.maxSpeed );
    %     % max running speed of pups varies quite a lot so not surenormalising by max speed = 40cm/s makes sense. Instead normalise by max speed for current trial?
    %     % speedMaps{i,1}(:,4) = (speedMaps{i,1}(:,2) - b(1) ) ./ (b(2) * max(speedMaps{i,1}(:,1),[],'omitnan'));
    %     % this can yield negative values if rate is lower than y intercept of regression - maybe set those to 0?
    % end

% maxSpeed          = round(prctile(obj.posData.speed{trialInd},99)); % Hasselmo
% [speedHist,bins,ind] = histcounts(obj.posData.speed{trialInd},linspace(0,obj.mapParams.speed.maxSpeed,ceil(obj.mapParams.speed.maxSpeed/obj.mapParams.speed.binSizeSpeed)+1));
% [speedHist,~,ind] = histcounts(obj.posData.speed{trialInd},linspace(0,maxSpeed,ceil(maxSpeed/obj.mapParams.speed.binSizeSpeed)+1));

% validBins         = speedHist > obj.mapParams.speed.minBinProp * sum(speedHist); % in Kropf et al. it's >0.5% of all bins
%
% occCounts         = accumarray(ind(ind~=0), 1, [length(speedHist) 1])./ obj.trialMetaData(trialInd).posFs;

% speedBins         = obj.mapParams.speed.binSizeSpeed/2:obj.mapParams.speed.binSizeSpeed:obj.mapParams.speed.maxSpeed-obj.mapParams.speed.binSizeSpeed/2;


    % spikeCounts             = accumarray(ind(ind~=0),instSpikeCount(ind~=0)',[length(speedHist) 1]);

 % spikeCounts(~validBins) = NaN;
    % maybe non-valid bins should be excluded from smoothing? 
    % speedMaps{i,1}(:,2)          = imfilter(spikeCounts./occCounts,fspecial('gaussian',[ceil(obj.mapParams.speed.smKernelLength/obj.mapParams.speed.binSizeSpeed) 1],3/obj.mapParams.speed.binSizeSpeed));
    % kernelSD                         = obj.mapParams.speed.smKernelLength/obj.mapParams.speed.binSizeSpeed; 
    % kernel                           = fspecial('gaussian',[round(kernelSD*4)+1 1],kernelSD);
    % speedMaps{i,1}(:,2)              = scanpix.helpers.smoothWithNaNs(spikeCounts./occCounts,kernel);

        % confInt                          = cell2mat( cellfun(@(x) CI(x,obj.mapParams.speed.confInt), accumarray(ind(ind~=0),instSpikeCount(ind~=0)',[length(speedHist) 1],@(x) {x}),'uni',0) ) .* obj.trialMetaData(trialInd).posFs;

         % speedMaps{i,1}(~validBins,[1 3]) = NaN;

             % speedMaps{i,2} = corr(obj.posData.speed{trialInd}(~isnan(obj.posData.speed{trialInd})),instFRate(~isnan(obj.posData.speed{trialInd})));
