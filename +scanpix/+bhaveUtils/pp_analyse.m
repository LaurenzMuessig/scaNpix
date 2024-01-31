function Tout = pp_analyse(dataObj)


minSpeed     = 2.5;
feederRadius = 50;
minTimeBetweenCrossings = 0.5;

% Set up results table %
scoreDum     = nan(1,2);
varList =   {
    'rat',         nan; ...
    'age',         nan; ...
    'exp_group',   nan; ...
    'dataset',     'string'; ...
    'trialID', cell(size(scoreDum)); ....

    'pathLength',   scoreDum; ...
    'meanSpeed',    scoreDum; ...
    'prcntRunning', scoreDum; ...

    'rewCentre',   cell(size(scoreDum)); ...
    'rewRadius',   scoreDum; ...
    
    'nZoneCross',    scoreDum; ...
    'zoneCrossings', cell(size(scoreDum)); 
    'nRewTrigg',      scoreDum; ...
    'rewTriggered', cell(size(scoreDum));

    'feederToRewZoneInd', cell(size(scoreDum)); 

    };

varList = varList';
Tout = cell2table( varList(2,:) );
Tout.Properties.VariableNames = varList(1,:);
% Tout.Properties.UserData = prms;


for i = 1:length(dataObj.trialNames)

    Tout.rat        = sscanf(dataObj.trialMetaData(i).animal,'%*c%d');
    Tout.age        = dataObj.trialMetaData(i).age;
    Tout.exp_group  = dataObj.trialMetaData(i).group;
    Tout.dataset    = dataObj.dataSetName;
    Tout.trialID{1,i} = dataObj.trialNames(i);
    % path data
    Tout.pathLength(i)   = sum(sqrt( diff(dataObj.posData.XY{i}(dataObj.posData.speed{i}>minSpeed,1)).^2 + diff(dataObj.posData.XY{i}(dataObj.posData.speed{i}>minSpeed,2)).^2 ),'omitnan') ./ dataObj.trialMetaData(i).ppm; % in m
    Tout.meanSpeed(i)    = mean(dataObj.posData.speed{i},'omitnan');
    Tout.prcntRunning(i) = sum(dataObj.posData.speed{i}>minSpeed,'omitnan') / length(dataObj.posData.speed{i});
    % scale reward locations (like position data)
    if dataObj.trialMetaData(i).PosIsScaled
        rewCentre          = dataObj.trialMetaData(i).rewardZoneCentre .* (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org);
        radius             = dataObj.trialMetaData(i).rewardZoneRadius * (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org); 
        feederPos          = [dataObj.trialMetaData(i).feederCoords(1:2:end)', dataObj.trialMetaData(i).feederCoords(2:2:end)'] .* (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org);
        feederRadiusScaled = feederRadius * (dataObj.trialMetaData(i).ppm/dataObj.trialMetaData(i).ppm_org); 
    end

    if dataObj.trialMetaData(i).PosIsFitToEnv{1,1}
        rewCentre = rewCentre - [dataObj.trialMetaData(i).PosIsFitToEnv{1,2}(1), dataObj.trialMetaData(i).PosIsFitToEnv{1,2}(2)];
        feederPos = feederPos - [dataObj.trialMetaData(i).PosIsFitToEnv{1,2}(1), dataObj.trialMetaData(i).PosIsFitToEnv{1,2}(2)];
    end
    % keep a record of the scaled reward zone properties
    Tout.rewCentre{i} = rewCentre;
    Tout.rewRadius(i) = radius;

    % dwell in reward zone
    inRewInd = (dataObj.posData.XY{i}(:,1) - rewCentre(1)).^2 + (dataObj.posData.XY{i}(:,2) - rewCentre(2)).^2 < radius.^2;
    % n Zone crossings
    zoneCrossings = [find(diff([0;inRewInd;0])==1),find(diff([0;inRewInd;0])==-1)];
    % join crossings that are too close in time 
    while any(zoneCrossings(2:end,1) - zoneCrossings(1:end-1,2)  < minTimeBetweenCrossings * dataObj.params('posFs') )
        idx = find(zoneCrossings(2:end,1) - zoneCrossings(1:end-1,2)  < minTimeBetweenCrossings * dataObj.params('posFs'),1,'first');
        zoneCrossings(idx,2) = zoneCrossings(idx+1,2);
        zoneCrossings(idx+1,:) = [];
    end

    Tout.nZoneCross(i)      = size(zoneCrossings,1);
    Tout.zoneCrossings{1,i} = zoneCrossings;
    % n rewards
    rewTriggered            = find(diff([0;double(dataObj.bhaveData.data{i});0])>0);
  
    % feeder visits
    atFeeder = false(length(dataObj.posData.XY{i}),1);
    nFeeder  = zeros(length(dataObj.posData.XY{i}),1);
    for j = 1:size(feederPos,1)    
        atFeeder = atFeeder | (dataObj.posData.XY{i}(:,1) - feederPos(j,1)).^2 + (dataObj.posData.XY{i}(:,2) - feederPos(j,2)).^2 < feederRadiusScaled.^2;
        nFeeder((dataObj.posData.XY{i}(:,1) - feederPos(j,1)).^2 + (dataObj.posData.XY{i}(:,2) - feederPos(j,2)).^2 < feederRadiusScaled.^2) = j;
    end
    % now get those zone crossings that triggered reward delivery (first
    % pilot only)
    zoneCrossingsRewTrig = nan(length(rewTriggered),2);
    c = 1;
    remInd = [];
    for j = 1:length(rewTriggered)
        [diffVal,tmpInd] = min(abs(zoneCrossings(:,1)-rewTriggered(j))); %
        if ~isempty(tmpInd) && diffVal < dataObj.params('posFs')
            zoneCrossingsRewTrig(c,:) = zoneCrossings(tmpInd,:);
            c = c+1;
        elseif diffVal > dataObj.params('posFs')
            % in case a manual reward was triggered by accident we want to
            % remove those from the list - happened at least once
            remInd = [remInd;j];
        end
    end
    rewTriggered(remInd) = [];
    zoneCrossingsRewTrig(isnan(zoneCrossingsRewTrig(:,1)),:) = [];

    Tout.nRewTrigg(i)       = length(rewTriggered);
    Tout.rewTriggered{1,i}  = rewTriggered;
    
    % 
    atFeederInd        = find(atFeeder);
    feederToRewZoneInd = false(length(dataObj.posData.XY{i}),1);
    for j = 1:size(zoneCrossingsRewTrig,1)
        tmpInd        = find(atFeederInd<zoneCrossingsRewTrig(j,1),1,'last');
        lastFeederInd = atFeederInd(tmpInd);
        feederToRewZoneInd(lastFeederInd:zoneCrossingsRewTrig(j,1)) = true;
    end
    Tout.feederToRewZoneInd{i} = feederToRewZoneInd;
end

end

