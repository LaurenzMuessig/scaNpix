function [coords, coords_binned, radius] = getBarrCoords(obj, trInd, barrType)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trInd (1,:) {mustBeNumericOrLogical}
    barrType (1,:) {mustBeMember(barrType,{'straight','circ'})} = 'straight';
end


%%
objFieldStrs = getValidObjStrings(obj.trialMetaData(trInd));

%%

if ~isempty(objFieldStrs)

    % fetch data - format depends on barrier type
    radius = nan(length(objFieldStrs),2);
    if strcmp(barrType,'straight')
        coords   = nan(length(objFieldStrs),8);
        circFlag = false;
    else
        coords   = nan(length(objFieldStrs),2);
        circFlag = true;
    end
    %
    for i = 1:length(objFieldStrs)
        if ~circFlag
            coords(i,:) = [obj.trialMetaData(trInd).(objFieldStrs{i})(1:2:end) obj.trialMetaData(trInd).(objFieldStrs{i})(2:2:end)];
        else
            coords(i,:) = [obj.trialMetaData(trInd).(objFieldStrs{i})(1) obj.trialMetaData(trInd).(objFieldStrs{i})(2)];
            radius(i,1) = obj.trialMetaData(trInd).(objFieldStrs{i})(3);
        end
    end
    % convert to frame of position data (scaling to common ppm and fit to environment)
    % coords are [x1..xn y1..yn] per row, so convert as points and put back into same format
    nPts                  = size(coords,2)/2;
    X                     = coords(:,1:nPts);
    Y                     = coords(:,nPts+1:end);
    [xyConv, lenScale]    = scanpix.maps.rawToFitCoords(obj, trInd, [X(:) Y(:)]);
    coords                = [reshape(xyConv(:,1),size(X)) reshape(xyConv(:,2),size(Y))];
    radius                = radius .* mean(lenScale); % fit to env. scaling can differ slightly between x/y
    % also generate a binned version of the corrdinates (bin size of rate maps)
    % need to replicate the frame used in scanpix.maps.makeRateMaps: if position data are neither fit to the environment
    % nor scaled to the camera window, the minimum (speed filtered) position is subtracted before binning
    if ~obj.trialMetaData(trInd).PosIsFitToEnv{1} && ~(isKey(obj.params,'scalePos2CamWin') && obj.params('scalePos2CamWin'))
        positions = obj.posData.XY{trInd};
        if obj.mapParams.rate.speedFilterFlagRMaps
            speedFilter              = obj.posData.speed{trInd} <= obj.mapParams.rate.speedFilterLimitLow | obj.posData.speed{trInd} > obj.mapParams.rate.speedFilterLimitHigh;
            positions(speedFilter,:) = NaN;
        end
        minXY     = min(positions,[],1,'omitnan');
        coordsMap = [coords(:,1:nPts)-minXY(1) coords(:,nPts+1:end)-minXY(2)];
    else
        coordsMap = coords;
    end
    % rate maps bin with ceil(XY/binSizePix), i.e. bin k spans ((k-1)*binSizePix, k*binSizePix], whereas vector maps
    % treat integers as bin centres, so shift by 1/2 bin
    binSizePix    = floor( obj.trialMetaData(trInd).ppm/100 * obj.mapParams.rate.binSizeSpat );
    coords_binned = coordsMap ./ binSizePix + 0.5;
    radius(:,2)   = radius(:,1) ./ binSizePix;
else
    [coords,coords_binned, radius] = deal([]);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [objFieldStrs] = getValidObjStrings(trialMetaStruct)

fNames       = fieldnames(trialMetaStruct);
ind          = ~cellfun('isempty',regexp(fNames,'objectPos(\d|)'));
objFieldStrs = fNames(ind);
%
allEmptyInd  = structfun(@isempty,trialMetaStruct);
emptyFields  = fNames(ind & allEmptyInd);
objFieldStrs = objFieldStrs(~ismember(objFieldStrs,emptyFields));

end