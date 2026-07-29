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
            radius(i)   = NaN;
        else
            coords(i,:) = [obj.trialMetaData(trInd).(objFieldStrs{i})(1) obj.trialMetaData(trInd).(objFieldStrs{i})(2)];
            radius(i)   = obj.trialMetaData(trInd).(objFieldStrs{i})(3);
        end
    end
    % add scaling factor in case data is scaled to common ppm
    if obj.trialMetaData(trInd).PosIsScaled
        scaleFact = obj.trialMetaData(trInd).ppm / obj.trialMetaData(trInd).ppm_org;
    else
        scaleFact = 1;
    end
    %
    coords = coords .* scaleFact;
    radius = radius .* scaleFact;
   
    % in case pos is fitted to visited environment or embedded in camera window need to adjust coordinates further
    if obj.trialMetaData(trInd).PosIsFitToEnv{1}
        coords = [coords(:,1:size(coords,2)/2) - obj.trialMetaData(trInd).PosIsFitToEnv{2}(1) coords(:,size(coords,2)/2+1:end) - obj.trialMetaData(trInd).PosIsFitToEnv{2}(2)];
    end
    % also generate a binned version of the corrdinates (bin size of rate maps)
    binSizePix    = floor( obj.trialMetaData(trInd).ppm/100 * obj.mapParams.rate.binSizeSpat );
    coords_binned = coords ./ binSizePix;
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