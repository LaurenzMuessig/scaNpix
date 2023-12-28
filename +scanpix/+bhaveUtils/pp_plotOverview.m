function pp_plotOverview(dataObj,ResTable)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if nargin < 2
    ResTable = scanpix.bhaveUtils.pp_analyse(dataObj);
end

% open figure
axArray = scanpix.plot.multPlot([length(dataObj.trialNames) 3]);
% make pos maps if necessary
if isempty(dataObj.maps.pos{1})
    dataObj.addMaps('pos');
end
% plot
for i = 1:length(dataObj.trialNames)
    
    plot(axArray{i,1},dataObj.posData.XY{i}(:,1),dataObj.posData.XY{i}(:,2),'k-');
    hold(axArray{i,1},'on');
    plot(axArray{i,1},sin(linspace(0,2*pi,50)).*ResTable.rewRadius(1,i)+ResTable.rewCentre{i}(1),cos(linspace(0,2*pi,50)).*ResTable.rewRadius(1,i)+ResTable.rewCentre{i}(2),'r-');
    hold(axArray{i,1},'off');
    %
    scanpix.plot.plotRateMap(dataObj.maps.pos{i}{1},axArray{i,2});
    hold(axArray{i,2},'on');
    rewCentreBinned = ResTable.rewCentre{i} ./ floor( dataObj.trialMetaData(i).ppm/100 * dataObj.mapParams.rate.binSizeSpat ); 
    radiusBinned    = ResTable.rewRadius(1,i) ./ floor( dataObj.trialMetaData(i).ppm/100 * dataObj.mapParams.rate.binSizeSpat ); 
    plot(axArray{i,2},sin(linspace(0,2*pi,50)).*radiusBinned+rewCentreBinned(1),cos(linspace(0,2*pi,50)).*radiusBinned+rewCentreBinned(2),'r-');
    hold(axArray{i,2},'off');
    %
    tmp = dataObj.posData.XY{i};
    tmp(~ResTable.feederToRewZoneInd{i},:) = NaN;
    plot(axArray{i,3},tmp(:,1),tmp(:,2),'k-');

end






end

