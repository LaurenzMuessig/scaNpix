function pp_plotOverview(dataObj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


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
    rewCentre = ((dataObj.trialMetaData(i).rewardZoneCentre.* (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org))-min(dataObj.posData.XYraw{i},[],1)); 
    radius    = (dataObj.trialMetaData(i).rewardZoneRadius * (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org)); 
    plot(axArray{i,1},sin(linspace(0,2*pi,50)).*radius+rewCentre(1),cos(linspace(0,2*pi,50)).*radius+rewCentre(2),'r-');
    hold(axArray{i,1},'off');
    %
    scanpix.plot.plotRateMap(dataObj.maps.pos{i}{1},axArray{i,2});
    hold(axArray{i,2},'on');
    rewCentreBinned = rewCentre ./ floor( dataObj.trialMetaData(i).ppm/100 * dataObj.mapParams.rate.binSizeSpat ); 
    radiusBinned    = radius ./ floor( dataObj.trialMetaData(i).ppm/100 * dataObj.mapParams.rate.binSizeSpat ); 
    plot(axArray{i,2},sin(linspace(0,2*pi,50)).*radiusBinned+rewCentreBinned(1),cos(linspace(0,2*pi,50)).*radiusBinned+rewCentreBinned(2),'r-');
    hold(axArray{i,2},'off');

end






end

