function [mapSeries, timeInt] = makeMapTimeSeries(obj,timeInt,trialInd,options)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    timeInt {mustBeNumeric}
    trialInd (1,1) {mustBeNumeric} = 1;
    options.type (1,:) {mustBeMember(options.type,{'rate','dir'})} = 'rate';
    options.rep (1,1) {mustBeNumericOrLogical} = true;
    options.cellind {mustBeNumericOrLogical} = true(size(obj.cell_ID,1),1);
    options.prms struct;
end

%
if isempty(options.prms)
    prms = obj.mapParams.(options.type);
else
    prms = options.prms;
end
%
if isKey(obj.params,'InterpPos2PosFs') && obj.params('InterpPos2PosFs')
    sampleTimes = [];
    prms.posFs  = obj.trialMetaData(trialInd).log.InterpPosFs;
else
    sampleTimes = obj.spikeData.sampleT{trialInd};
end
%
if strcmp(options.type,'rate') && ~isempty( obj.trialMetaData(trialInd).envSize ) 
    prms.envSize = obj.trialMetaData(trialInd).envSize ./ 100 .* obj.trialMetaData(trialInd).ppm; % in pixels
end


%%
if options.rep
    % expand time series in case time interval should be repeated
    timeInt = [ ( timeInt(1):timeInt(2)-timeInt(1):obj.trialMetaData(trialInd).duration-timeInt(2) )', ( timeInt(2):timeInt(2)-timeInt(1):obj.trialMetaData(trialInd).duration )' ];
end
%
mapSeries = cell(1,size(timeInt,1));
for i = 1:size(timeInt,1)
    keepInd                  = false(size(obj.posData.XY{trialInd},1),1);
    startInd                 = max([1,ceil(prms.posFs  * timeInt(i,1))]);
    endInd                   = ceil(prms.posFs * timeInt(i,2));
    keepInd(startInd:endInd) = true;

    switch options.type
        case 'rate'
            tempPos             = obj.posData.XY{trialInd};
            tempPos(~keepInd,:) = NaN;
            mapSeries{1,i}      = scanpix.maps.makeRateMaps(obj.spikeData.spk_Times{trialInd}(options.cellind),tempPos,sampleTimes,obj.trialMetaData(trialInd).ppm,obj.posData.speed{trialInd},prms);
        case 'dir'
            tempDir             = obj.posData.dir{trialInd};
            tempDir(~keepInd)   = NaN;
            mapSeries{1,i}      = scanpix.maps.makeDirMaps(obj.spikeData.spk_Times{trialInd}(options.cellind),tempDir,sampleTimes,obj.posData.speed{trialInd},prms);
        case 'speed'
            % to do
    end
end

end