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
if ~isempty(obj.trialMetaData(trInd).objectPos)

    % fetch data - format depends on barrier type
    if strcmp(barrType,'straight')
        coords = [obj.trialMetaData(trInd).objectPos(1:2:end)' obj.trialMetaData(trInd).objectPos(2:2:end)'];
        radius = NaN;
    else
        coords = [obj.trialMetaData(trInd).objectPos(1) obj.trialMetaData(trInd).objectPos(2)];
        radius = obj.trialMetaData(trInd).objectPos(3);
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
        coords = coords - obj.trialMetaData(trInd).PosIsFitToEnv{2};
    
    elseif obj.params('scalePos2CamWin')
        xMin = obj.trialMetaData(trInd).xmin * scaleFact;
        yMin = obj.trialMetaData(trInd).ymin * scaleFact;
        %
        coords = coords - [xMin yMin];
    end
    % also generate a binned version of the corrdinates (bin size of rate maps)
    binSizePix    = floor( obj.trialMetaData(trInd).ppm/100 * obj.mapParams.rate.binSizeSpat );
    coords_binned = coords ./ binSizePix;
    radius(1,2)   = radius ./ binSizePix;
else
    [coords,coords_binned, radius] = deal([]);
end


end