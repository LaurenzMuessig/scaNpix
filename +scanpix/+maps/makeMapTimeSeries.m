function [mapSeries, timeInt] = makeMapTimeSeries(obj,timeInt,trialInd,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% defaults
mapType        = 'rate';
repeatInterval = true;
cellInd        = true(size(obj.cell_ID,1),1);

%%
p = inputParser;
addOptional(p,'cellind', cellInd,        ( @(x) islogical(x) || isnumeric(x) ) );
addParameter(p,'type',   mapType,        ( @(x) ischar(x) || isstring(x) )     );
addParameter(p,'rep',    repeatInterval, @islogical                            );
%
parse(p,varargin{:});

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
            mapSeries{1,i}      = scanpix.maps.makeRateMaps(obj.spikeData.spk_Times{trialInd}(p.Results.cellind),tempPos,obj.spikeData.sampleT{trialInd},obj.trialMetaData.ppm,tempSpeed,obj.mapParams.rate);
        case 'dir'
            tempDir             = obj.posData.dir{trialInd};
            tempDir(~keepInd,:) = NaN;
            mapSeries{1,i}      = scanpix.maps.makeDirMaps(obj.spikeData.spk_Times{trialInd}(p.Results.cellind),tempDir,obj.spikeData.sampleT{trialInd},tempSpeed,obj.mapParams.dir);
        case 'speed'
            % to do
    end
end

end