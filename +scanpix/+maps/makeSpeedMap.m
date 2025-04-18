function speedMaps = makeSpeedMap(obj, trialInd, cellInd)
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
end

%%
spkTimes = obj.spikeData.spk_Times{trialInd}(cellInd);

%%
% anon. fcn for 95% conf interval calculation - assuming normal dist which
% is prob. not the best way
CI = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x(:)))-1); %%

%% make speed map

% bin running speeds
[instSpeed,~,ind] = histcounts(obj.posData.speed{trialInd},0:obj.mapParams.speed.binSizeSpeed:obj.mapParams.speed.maxSpeed);
validBins         = instSpeed > obj.mapParams.speed.minBinProp * obj.trialMetaData(trialInd).duration * obj.trialMetaData(trialInd).posFs; % in Kropf et al. it's >0.5% of all bins
%
occCounts         = accumarray(ind(ind~=0), 1, [size(instSpeed,2) 1])./ obj.trialMetaData(trialInd).posFs;
%
speedMaps         = cell(length(obj.spikeData.spk_Times{trialInd}),3);

if obj.mapParams.speed.showWaitBar; hWait = waitbar(0); end

for i = 1:length(spkTimes)
    % inst firing rate
    instFRate = histcounts(ceil(spkTimes{i} * obj.trialMetaData(trialInd).posFs),0:length(obj.posData.speed{trialInd})); %.* prms.posFs;
    % smooth
    kernel    = ones(1,ceil(obj.mapParams.speed.smKernelLength * obj.trialMetaData(trialInd).posFs)) ./ (obj.mapParams.speed.smKernelLength * obj.trialMetaData(trialInd).posFs); % Kropf uses kernel of 250ms
%     h = fspecial('average',[ceil(prms.smKernelLength * prms.posFs) 1],1/prms.posFs)
    instFRate = imfilter(instFRate,kernel,'replicate');
    
    % mean rate / speed bin 
    speedMaps{i,1}(:,1) = obj.mapParams.speed.binSizeSpeed/2:obj.mapParams.speed.binSizeSpeed:obj.mapParams.speed.maxSpeed-obj.mapParams.speed.binSizeSpeed/2;
    spikeCounts         = accumarray(ind(ind~=0),instFRate(ind~=0)',[length(instSpeed) 1]);
    speedMaps{i,1}(:,2) = spikeCounts./occCounts;
%     speedMaps{i,1}(:,2) = accumarray(ind(ind~=0),instFRate(ind~=0)',[length(instSpeed) 1],@mean);
    % 95 CI for means (a bit dense)
%     speedMaps{i,1}(:,3:4) = cell2mat( cellfun(@(x) CI(x,prms.confInt), accumarray(ind(ind~=0),instFRate(ind~=0)',[length(instSpeed) 1],@(x) {x}),'uni',0) );
    speedMaps{i,1}(:,3:4) = cell2mat( cellfun(@(x) CI(x,obj.mapParams.speed.confInt), accumarray(ind(ind~=0),instFRate(ind~=0)',[length(instSpeed) 1],@(x) {x}),'uni',0) ) .* obj.trialMetaData(trialInd).posFs;
    speedMaps{i,1}        = speedMaps{i,1}(validBins,:);
    % add Kropf et al. normalisation
%     if prms.normaliseFR
%         b =  regress(instFRate',[ind ones(length(ind), 1)]);
% %         yInt = speedMaps{i}(1,2) - b(1)*speedMaps{i}(1,1);
%         tmpMap         = accumarray(ind(ind~=0),(instFRate(ind~=0)'-b(2))./(b(1)*prms.maxSpeed),[length(instSpeed) 1],@mean);
%         speedMaps{i,2} = tmpMap(validBins);
%         speedMaps{i,3} = [b(1) b(2)];
%     end
    if obj.mapParams.speed.showWaitBar; waitbar(i/length(spkTimes),hWait,sprintf('Making those Speed Maps... %i/%i done.',i,length(spkTimes))); end

end

if obj.mapParams.speed.showWaitBar; close(hWait); end

end

