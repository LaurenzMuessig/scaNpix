function [mapSeries, timeInt] = makeMapTimeSeries(obj,timeInt,trialInd,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% defaults
mapType        = 'rate';
repeatInterval = true;
cellInd        = true(size(obj.cell_ID,1),1);
mapParams      = scanpix.maps.defaultParamsRateMaps;

%%
p = inputParser;
addOptional(p, 'cellind', cellInd,        ( @(x) islogical(x) || isnumeric(x) ) );
addParameter(p,'type',    mapType,        ( @(x) mustBeMember(x,{'rate','dir'})));
addParameter(p,'rep',     repeatInterval, @islogical                            );
addOptional(p, 'prms',    mapParams,      ( @(x) isstruct(x) || isempty(x) )    );
%
parse(p,varargin{:});
%
if isempty(p.Results.prms)
    prms = obj.mapParams.(p.Results.type);
else
    prms = p.Results.prms;
end
%
if isKey(obj.params,'InterpPos2PosFs') && obj.params('InterpPos2PosFs')
    sampleTimes = [];
    prms.posFs  = obj.trialMetaData(trialInd).log.InterpPosFs;
else
    sampleTimes = obj.spikeData.sampleT{trialInd};
end
%
if strcmp(p.Results.type,'rate') && ~isempty( obj.trialMetaData(trialInd).envSize ) 
    prms.envSize = obj.trialMetaData(trialInd).envSize ./ 100 .* obj.trialMetaData(trialInd).ppm; % in pixels
end


%%
if p.Results.rep
    % expand time series in case time interval should be repeated
    timeInt = [ ( timeInt(1):timeInt(2)-timeInt(1):obj.trialMetaData(trialInd).duration-timeInt(2) )', ( timeInt(2):timeInt(2)-timeInt(1):obj.trialMetaData(trialInd).duration )' ];
end
%
mapSeries = cell(1,size(timeInt,1));
for i = 1:size(timeInt,1)
    keepInd                  = false(size(obj.posData.XY{trialInd},1),1);
    startInd                 = max([1,obj.params('posFs') * timeInt(i,1)]);
    endInd                   = obj.params('posFs') * timeInt(i,2);
    keepInd(startInd:endInd) = true;
    tempSpeed                = obj.posData.speed{trialInd};
    tempSpeed(~keepInd,:)    = NaN;

    switch p.Results.type
        case 'rate'
            tempPos             = obj.posData.XY{trialInd};
            tempPos(~keepInd,:) = NaN;
            mapSeries{1,i}      = scanpix.maps.makeRateMaps(obj.spikeData.spk_Times{trialInd}(p.Results.cellind),tempPos,sampleTimes,obj.trialMetaData(trialInd).ppm,tempSpeed,prms);
        case 'dir'
            tempDir             = obj.posData.dir{trialInd};
            tempDir(~keepInd,:) = NaN;
            mapSeries{1,i}      = scanpix.maps.makeDirMaps(obj.spikeData.spk_Times{trialInd}(p.Results.cellind),tempDir,sampleTimes,tempSpeed,prms);
        case 'speed'
            % to do
    end
end

end