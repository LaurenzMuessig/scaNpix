function [] = pp_analyse(dataObj)
%UNTITLED4 Summary of this function goes here


for i = 1:length(dataObj.trialNames)
    
    rewCentre = ((dataObj.trialMetaData(i).rewardZoneCentre.* (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org))-min(dataObj.posData.XYraw{i},[],1)); 
    radius    = (dataObj.trialMetaData(i).rewardZoneRadius * (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org)); 
    % dwell in reward zone
    dwellInd = (dataObj.posData.XY{i}(:,1) - rewCentre(1)).^2 + (dataObj.posData.XY{i}(:,2) - rewCentre(2)).^2 <= radius.^2;
    

end



end

