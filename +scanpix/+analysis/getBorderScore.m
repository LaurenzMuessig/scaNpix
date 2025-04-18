function [borderScore, mainWall] = getBorderScore(maps,binSize,options)
% Calculate border-score, as defined by Solstad et al, 2008. This is an
% updated version from TW's original function ('borderScoreForScaledRateMaps'). 
% This version doesn't assume the rate map is a square array, but will
% define the border of the map by the sampling of the rate map. The box
% boundaries for each wall are then set as the average extent of the map 
% in each direction. 
% Note: Any full row/column of the map that is only nan's will be trimmed
% off. For a perfectly sampled square map, it will return same border
% score as original function
%
% Usage:
%
%       [BS, mainWall] = getBorderScore(maps, binSizeInCmSq);
%       [BS, mainWall] = getBorderScore(maps, binSizeInCmSq, optionalInputStruct );
%       [BS, mainWall] = getBorderScore(maps, binSizeInCmSq, 'inputName', inputVal, .. etc .. );
%
% Inputs:   maps          - rate map / cell array of maps
%           binSize       - bin size of rate map in cm^2
%
%
% Options:
%
%   'rateThr',          0.2,    - Threshold for defining fields, proportion of peak rate.
%   'sizeThr',          200,    - Threshold for minimum field size to enter into the analysis.
%   'maxBinDistBorder', 4       - maximum wall distance (from max extent) for a bin to be counted as part of the border of box
%   'debugModeON',      0       - makes a debug plot; rate map and field map with Cm bins marked 
%
%
% Outputs:  borderScore   - border score for rate map
%           mainWall      - index of wall with largest field coverage (Cm);
%                           1=N, 2=E, 3=S, 4=W
%
% TW, LM 2019

%% 
%%
arguments
    maps 
    binSize (1,1) {mustBeNumeric}
    options.rateThr (1,1) {mustBeNumeric} = 0.2;
    options.sizeThr (1,1) {mustBeNumeric} = 200;
    options.maxBinDistBorder (1,1) {mustBeNumeric} = 4;
    options.debug (1,1) {mustBeNumericOrLogical} = false;
end
                                                                             %
%%
if ~iscell(maps)
    maps = {maps};
end

% trim of any margin that is purely NaN
nanIndC1 = find(sum(isnan(maps{1}),1) ~= size(maps{1},1),1,'first');
nanIndC2 = find(sum(isnan(maps{1}),1) ~= size(maps{1},1),1,'last');
nanIndR1 = find(sum(isnan(maps{1}),2) ~= size(maps{1},2),1,'first');
nanIndR2 = find(sum(isnan(maps{1}),2) ~= size(maps{1},2),1,'last');
maps     = cellfun(@(x) x(nanIndR1:nanIndR2,nanIndC1:nanIndC2), maps, 'uni', 0);
%
templateMap = maps{1};
%% FIND BORDERS
% N.B. could also use 'bwboundaries'

[cols,rows] = meshgrid(1:size(templateMap,2),1:size(templateMap,1)); % map coord grid
colSub      = 1:size(templateMap,2); % index for all columns 
rowSub      = 1:size(templateMap,1); % index for all rows 
% bins for borders in all 4 directions
% North
indN                = accumarray(cols(:),templateMap(:),[size(templateMap,2) 1],@(x) find(~isnan(x),1,'first')); % this is the row index of northern map edge (orderd by column N)
filtBinOutInd       = indN > options.maxBinDistBorder;
indN(filtBinOutInd) = []; % remove bins that have too high wall distance
meanN               = floor( mean(indN) ); % average
indN                = [indN colSub(~filtBinOutInd)']; % final index into map for border bins
% now do the rest for other walls
% south
indS                = accumarray(cols(:),maps{1}(:),[size(maps{1},2) 1],@(x) find(~isnan(x),1,'last'));
filtBinOutInd       = indS < max(indS) - options.maxBinDistBorder;
indS(filtBinOutInd) = [];
meanS               = ceil(mean(indS));
indS                = [indS colSub(~filtBinOutInd)'];
% for E/W we need to transpose everything
rows                = rows';
templateMap         = templateMap';
indE                = accumarray(rows(:),templateMap(:),[size(templateMap,2) 1],@(x) find(~isnan(x),1,'first')); % this is the column index of Eastern map edge (orderd by row N)
filtBinOutInd       = indE > options.maxBinDistBorder;
indE(filtBinOutInd) = [];
meanE               = floor(mean(indE));
indE                = [rowSub(~filtBinOutInd)' indE];
% west
indW                = accumarray(rows(:),templateMap(:),[size(templateMap,2) 1],@(x) find(~isnan(x),1,'last'));
filtBinOutInd       = indW < max(indW) - options.maxBinDistBorder;
indW(filtBinOutInd) = [];
meanW               = ceil(mean(indW));
indW                = [rowSub(~filtBinOutInd)' indW];

templateMap = templateMap'; % transpose back to original format

%% Find Fields 
[borderScore, mainWall] = deal(nan(length(maps),1));
for i = 1:length(maps)
    
    % find all fields above thresh
    bwMap    = maps{i} >= ( max(maps{i}(:),[],'omitnan')*options.rateThr );
    fieldMap = bwlabel(bwMap,4);
    stats    = regionprops(fieldMap,'Area','PixelList');
    % only keep fields that are big enough
    keepInd  = [ stats(:).Area ] .* binSize > options.sizeThr ;
    % fail gracefully
    if ~any(keepInd); continue; end
    
    stats           = stats(keepInd); % remove too small fields
    % list of all field bins that pass criterion (rate + size)
    binListFields   = fliplr(vertcat( stats(:).PixelList )); % combined bin indices of all fields
    
    %% get Cm
    % calculate coverage on all 4 walls
    coverageOnFourWalls(1) = sum(ismember(binListFields,indN,'rows'));
    coverageOnFourWalls(2) = sum(ismember(binListFields,indE,'rows'));
    coverageOnFourWalls(3) = sum(ismember(binListFields,indS,'rows'));
    coverageOnFourWalls(4) = sum(ismember(binListFields,indW,'rows'));
    % calculate Cm
    [maxCov, mainWall(i)]  = max( coverageOnFourWalls, [], 'omitnan' );
    if ~mod(mainWall(i),2)
        Cm = maxCov / size(templateMap,1);
    else
        Cm = maxCov / size(templateMap,2);
    end
    
    %% get Dm
    % first we set the borders of box to be of size between average wall positions
    sizeX     = meanS - meanN + 1; % size box X
    sizeY     = meanW - meanE + 1; % size box Y
    hMapSizeX = ceil(sizeX/2);
    hMapSizeY = ceil(sizeY/2);
    % 4 cases (if not a square)
    if rem(sizeX,2) && ~rem(sizeY,2)
        [Y,X] = meshgrid([1:hMapSizeY fliplr(1:hMapSizeY)],[1:hMapSizeX fliplr(1:hMapSizeX-1)]);
    elseif ~rem(sizeX,2) && rem(sizeY,2)
        [Y,X] = meshgrid([1:hMapSizeY fliplr(1:hMapSizeY-1)],[1:hMapSizeX fliplr(1:hMapSizeX)]);
    elseif ~rem(sizeX,2) && ~rem(sizeY,2)
        [Y,X] = meshgrid([1:hMapSizeY fliplr(1:hMapSizeY)],[1:hMapSizeX fliplr(1:hMapSizeX)]);
    else
        [Y,X] = meshgrid([1:hMapSizeY fliplr(1:hMapSizeY-1)],[1:hMapSizeX fliplr(1:hMapSizeX-1)]);
    end
    dist2Wall = squeeze(min(  cat(3,X,Y)  ,[],  3  ));
    
    % now we need final distance map with all pixels outside of the bounded
    % distance map set to 1
    dist2wallFull                          = ones(size(templateMap));
    dist2wallFull(meanN:meanS,meanE:meanW) = dist2Wall; % insert distance map
    % Multiply this by the normalised FR map, and then get the mean distance to wall for the 'big enough' field pixels %
    % get firing rate of field pixels
    mapInd            = sub2ind(size(templateMap),binListFields(:,1),binListFields(:,2));
    rateInFields      = maps{i}(mapInd);
    % get wall distance of field pixels
    dist2WallInFields = dist2wallFull(mapInd);
    Dm                = transpose(dist2WallInFields) * (rateInFields / sum(rateInFields));
    % normnalise by half size of box
    Dm                = Dm / mean([hMapSizeY hMapSizeX]);
    
    %% border score
    % Calculate the border score, (Cm-Dm) / (Cm + Dm) %
    borderScore(i)    = (Cm-Dm) / (Cm + Dm);
    
    %% debug plot
    if options.debug
        figure;
        scanpix.plot.plotRateMap(maps{i},subplot(1,2,1));
        axis square
        fieldMap2PLot = ones(size(templateMap));
        fieldMap2PLot(isnan(templateMap)) = 0;
        fieldMap2PLot(mapInd) = 2;
        subplot(1,2,2)
        imagesc(fieldMap2PLot); colormap(gca,[1 1 1; 0 0 1; 1 0 0]);
        axis square
        hold on
        if mainWall == 1
            tempInd = ismember(binListFields,indN,'rows');
            scatter(binListFields(tempInd,2),binListFields(tempInd,1),'filled');
        elseif mainWall == 2
            tempInd = ismember(binListFields,indE,'rows');
            scatter(binListFields(tempInd,2),binListFields(tempInd,1),'filled');
        elseif mainWall == 3
            tempInd = ismember(binListFields,indS,'rows');
            scatter(binListFields(tempInd,2),binListFields(tempInd,1),'filled');
        else
            tempInd = ismember(binListFields,indW,'rows');
            scatter(binListFields(tempInd,2),binListFields(tempInd,1),'filled');
        end
        hold off
    end 
end

end

