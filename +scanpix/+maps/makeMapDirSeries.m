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
    dirRad        = mod(atan2(diff(obj.posData.XY{trialInd}(:,2)), diff(obj.posData.XY{trialInd}(:,1))),2*pi); %
    dirRad(end+1) = dirRad(end);
end

%%

% orgSt = obj.spikeData.spk_Times{trialInd};     

mapSeries = cell(1,length(dirShift));
% prevInd = false(length(dirRad),1);
for i = 1:length(dirShift)

    keepInd                      = any(cos(dirRad - dirShift(:,i)') > cos(dirRange/2),2);

    switch options.type
        case 'rate'
            mapSeries{1,i}       = scanpix.maps.makeRateMaps(obj, trialInd, 'addPosFilter',~keepInd,'cellInd', options.cellind );
            % mapSeries{1,i}       = scanpix.maps.makeRateMaps(obj, trialInd,'cellInd', options.cellind );
        case 'pos'
            [~,mapSeries{1,i},~] = scanpix.maps.makeRateMaps(obj, trialInd, 'mapType', 'pos', 'addPosFilter', ~keepInd );
        case 'speed'
            % to do
    end
    % obj.spikeData.spk_Times{trialInd} = orgSt;
end

end

    % keepInd = randperm(length(dirRad),sum(keepInd));
    % 
    % for j = 1:length(obj.spikeData.spk_Times{trialInd})
    % 
    %     st_filt{j} = scanpix.helpers.getFilteredSTimes(obj.spikeData.spk_Times{trialInd}{j},keepInd,obj.trialMetaData(trialInd).posFs); 
    % end
    % 
    % prevInd
    % obj.spikeData.spk_Times{trialInd} = st_filt;
    % numInd = 1:length(dirRad);
    % numInd = numInd(~prevInd);
    % tmp = zeros(length(dirRad),1);
    % ind = randperm(length(numInd),sum(keepInd));
    % tmp(numInd(ind)) = true;
    % keepInd = tmp;

    % prevInd = zeros(length(dirRad),1);
    % prevInd(numInd(ind)) = true;