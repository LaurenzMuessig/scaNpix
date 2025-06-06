function [mapSeries,dirShift] = makeMapDirSeries(obj,dirShift,dirRange,trialInd,options)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    dirShift {mustBeNumeric} = linspace(0,350,36) .* pi/180;
    dirRange (1,1) {mustBeNumeric} = 150 * pi/180; % use separate rows for miltiple dir centres
    trialInd (1,1) {mustBeNumeric} = 1;
    options.useHD (1,1) {mustBeNumericOrLogical} = true;
    options.type (1,:) {mustBeMember(options.type,{'rate','pos'})} = 'rate';
    options.cellind {mustBeNumericOrLogical} = true(size(obj.cell_ID,1),1);
end

%%
if options.useHD
    dirRad = obj.posData.direction{trialInd} .* pi/180;
else
    % movement direction
    dirRad        = mod(atan2(diff(obj.posData.XYraw{trialInd}(:,2)), diff(obj.posData.XYraw{trialInd}(:,1))),2*pi); %
    dirRad(end+1) = dirRad(end);
end

%%
mapSeries = cell(1,length(dirShift));
for i = 1:length(dirShift)

    keepInd                      = any(cos(dirRad - dirShift(:,i)') > cos(dirRange/2),2);

    switch options.type
        case 'rate'
            mapSeries{1,i}       = scanpix.maps.makeRateMaps(obj, trialInd, 'addPosFilter',~keepInd,'cellInd', options.cellind );
        case 'pos'
            [~,mapSeries{1,i},~] = scanpix.maps.makeRateMaps(obj, trialInd, 'mapType', 'pos', 'addPosFilter', ~keepInd );
        case 'speed'
            % to do
    end
end

end