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
    startInd                 = max([1,ceil( obj.trialMetaData(trialInd).posFs  * timeInt(i,1))]);
    endInd                   = ceil( obj.trialMetaData(trialInd).posFs * timeInt(i,2));
    keepInd(startInd:endInd) = true;

    switch options.type
        case 'rate'
            mapSeries{1,i}      = scanpix.maps.makeRateMaps(obj, trialInd, 'rate', ~keepInd, options.cellind );
        case 'dir'
            mapSeries{1,i}      = scanpix.maps.makeDirMaps(obj, trialInd, ~keepInd, options.cellind );
        case 'speed'
            % to do
    end
end

end