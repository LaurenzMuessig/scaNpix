function [borderScore, mainWall] = getBorderScore(map,binSize,varargin)
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
%       [BS, mainWall] = map_bordercell(map, binSizeInCmSq);
%       B = map_bordercell(map, binSizeInCmSq, optionalInputStruct );
%       B = map_bordercell(map, binSizeInCmSq, 'inputName', inputVal, .. etc .. );
%
% Inputs:   map           - rate map
%           binSizeInCmSq - bin size of rate map in cm^2
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

%% PARAMS
% Important variables %
prms.rateThr            = 0.2;       % Fraction of max rate to use as field threshold.
prms.sizeThr            = 200;       % This many cm-sq of contiguous bins.
prms.maxBinDistBorder   = 4;
prms.debugModeON        = 0; 

% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %

%% Pre Process map
% exit straight away if map is empty
if all( isnan(map(:)) )
    borderScore = nan;
    mainWall    = nan;
    return
end   
% trim of any margin that is purely NaN
nanIndC1 = find(sum(isnan(map),1) ~= size(map,1),1,'first');
nanIndC2 = find(sum(isnan(map),1) ~= size(map,1),1,'last');
nanIndR1 = find(sum(isnan(map),2) ~= size(map,2),1,'first');
nanIndR2 = find(sum(isnan(map),2) ~= size(map,2),1,'last');
map      = map(nanIndR1:nanIndR2,nanIndC1:nanIndC2);

%% FIND BORDERS
% N.B. could also use 'bwboundaries'
[cols,rows] = meshgrid(1:size(map,2),1:size(map,1)); % map coord grid
colSub      = 1:size(map,2); % index for all columns 
rowSub      = 1:size(map,1); % index for all rows 
% bins for borders in all 4 directions
% North
indN                = accumarray(cols(:),map(:),[size(map,2) 1],@(x) find(~isnan(x),1,'first')); % this is the row index of northern map edge (orderd by column N)
filtBinOutInd       = indN > prms.maxBinDistBorder;
indN(filtBinOutInd) = []; % remove bins that have too high wall distance
meanN               = floor( mean(indN) ); % average
indN                = [indN colSub(~filtBinOutInd)']; % final index into map for border bins
% now do the rest for other walls
% south
indS                = accumarray(cols(:),map(:),[size(map,2) 1],@(x) find(~isnan(x),1,'last'));
filtBinOutInd       = indS < max(indS) - prms.maxBinDistBorder;
indS(filtBinOutInd) = [];
meanS               = ceil(mean(indS));
indS                = [indS colSub(~filtBinOutInd)'];
% for E/W we need to transpose everything
rows                = rows';
map                 = map';
indE                = accumarray(rows(:),map(:),[size(map,2) 1],@(x) find(~isnan(x),1,'first')); % this is the column index of Eastern map edge (orderd by row N)
filtBinOutInd       = indE > prms.maxBinDistBorder;
indE(filtBinOutInd) = [];
meanE               = floor(mean(indE));
indE                = [rowSub(~filtBinOutInd)' indE];
% west
indW                = accumarray(rows(:),map(:),[size(map,2) 1],@(x) find(~isnan(x),1,'last'));
filtBinOutInd       = indW < max(indW) - prms.maxBinDistBorder;
indW(filtBinOutInd) = [];
meanW               = ceil(mean(indW));
indW                = [rowSub(~filtBinOutInd)' indW];

map = map'; % transpose back to original format

%% Find Fields 
% find all fields above thresh
bwMap    = map >= ( nanmax(map(:))*prms.rateThr );
fieldMap = bwlabel(bwMap,4);
stats    = regionprops(fieldMap,'Area','PixelList');
% only keep fields that are big enough
keepInd  = [ stats(:).Area ] .* binSize > prms.sizeThr ;
% fail gracefully
if ~any(keepInd)
    borderScore = nan;
    mainWall    = nan;
    return
end
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
[maxCov, mainWall]     = nanmax( coverageOnFourWalls );
if ~mod(mainWall,2)
    Cm = maxCov / size(map,1);
else
    Cm = maxCov / size(map,2);
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
dist2wallFull                          = ones(size(map));
dist2wallFull(meanN:meanS,meanE:meanW) = dist2Wall; % insert distance map
% Multiply this by the normalised FR map, and then get the mean distance to wall for the 'big enough' field pixels %
% get firing rate of field pixels
mapInd            = sub2ind(size(map),binListFields(:,1),binListFields(:,2));
rateInFields      = map(mapInd);
% get wall distance of field pixels
dist2WallInFields = dist2wallFull(mapInd);
Dm                = transpose(dist2WallInFields) * (rateInFields / sum(rateInFields));
% normnalise by half size of box
Dm                = Dm / mean([hMapSizeY hMapSizeX]);

%% border score
% Calculate the border score, (Cm-Dm) / (Cm + Dm) %
borderScore       = (Cm-Dm) / (Cm + Dm);

%% debug plot
if prms.debugModeON
    figure;
    subplot(1,2,1)
    gra_plotmap(map);
    fieldMap2PLot = ones(size(map));
    fieldMap2PLot(isnan(map)) = 0;
    fieldMap2PLot(mapInd) = 2;
    subplot(1,2,2)
    imagesc(fieldMap2PLot); colormap([1 1 1; 0 0 1; 1 0 0]);
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

