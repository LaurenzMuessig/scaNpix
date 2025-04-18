function [peakStats, peakMask] = fieldDetect(map,options)
%UNTITLED2 Summary of this function goes 
%   Detailed explanation goes here


%% TO DO

%%
arguments
    map {mustBeNumeric} 
    options.binMap (1,1) {mustBeNumericOrLogical} = true;
    options.nBinSteps (1,1) {mustBeNumeric} = 11; 
    options.binEdges (1,2) {mustBeNumeric} = [0 max(map(:),[],'omitnan')];
    options.thrMode (1,:) {mustBeMember(options.thrMode,{'abs','rel','none'})} = 'none';
    options.thr {mustBeScalarOrEmpty} = []; 
    options.minWSFieldSz (1,1) {mustBeNumeric} = 16; 
    options.minPeakSz (1,1) {mustBeNumeric} = 8; 
    options.debugOn (1,1) {mustBeNumericOrLogical} = true;
end


%%
if options.binMap
    tmpMap = scanpix.maps.binAnyRMap(map,'nsteps',options.nBinSteps,'cmapEdge',options.binEdges);
else
    tmpMap = map;
end

% regionfill(map,isnan(map));


% mask unvisited positions
tmpMap(isnan(map)) = -Inf;
% and if desired, background pixels
switch options.thrMode
    case 'abs'
        tmpMap(tmpMap < options.thr) = -Inf;
    case 'rel'
        tmpMap(tmpMap < max(tmpMap(:),[],'omitnan') * options.thr) = -Inf;
end

% segment the map using watershed
fieldsLabel = watershed(-tmpMap);
% special case if only a single field is present - watershed returns all 1's in that cxase
if all(fieldsLabel(:))
    fieldsLabel = tmpMap ~= -Inf;
end

%% merge fields that are too small into the closest larger one
fieldMask    = fieldsLabel ~= 0;
tmpStats     = regionprops(fieldMask,fieldsLabel,'centroid','PixelList','MaxIntensity');
% sort by size
[sz,sortInd] = sort( arrayfun(@(x) size(x.PixelList,1), tmpStats) );
tmpStats     = tmpStats(sortInd);
% index of fields with too small size
tooSmallInd  = find(sz < options.minWSFieldSz)';

% all field boundary pixels
fieldBorderPix = arrayfun(@(x) bwdist(fieldsLabel==x.MaxIntensity) < 2 & ~(fieldsLabel==x.MaxIntensity),tmpStats, 'UniformOutput',0);

while ~isempty(tooSmallInd)
    structInd         = 1:length(tmpStats);
    % all field centroid distances
    centroids         = reshape([tmpStats.Centroid],2,[])';
    dists             = squareform(pdist(centroids));
    dists(dists == 0) = NaN;
    % closest field
    [~,minInd]        = min(dists(tooSmallInd(1),:),[],'omitnan');
    % 
    otherFieldsInd    = ~ismember(structInd,[tooSmallInd(1);minInd]);
    mergeBorderPix    = fieldBorderPix{tooSmallInd(1)} & fieldBorderPix{minInd} & ~any(cat(3,fieldBorderPix{otherFieldsInd}),3);
    % update fields label
    fieldsLabel(mergeBorderPix)                                       = tmpStats(minInd).MaxIntensity;
    fieldsLabel(fieldsLabel == tmpStats(tooSmallInd(1)).MaxIntensity) = tmpStats(minInd).MaxIntensity;
    %
    % update structure after merging
    tmpStats(minInd).PixelList     = [tmpStats(minInd).PixelList; tmpStats(tooSmallInd(1)).PixelList];
    tmpStats(minInd).Centroid      = mean([tmpStats(minInd).Centroid; tmpStats(tooSmallInd(1)).Centroid],1);
    tmpStats(tooSmallInd(1))       = [];
    % update border mask after merging
    fieldBorderPix{minInd}         = (fieldBorderPix{tooSmallInd(1)} | fieldBorderPix{minInd}) & ~mergeBorderPix; 
    fieldBorderPix(tooSmallInd(1)) = [];
    % check for more fields < size thresh
    tooSmallInd                    = find(arrayfun(@(x) size(x.PixelList,1), tmpStats) < options.minWSFieldSz)';
end
 
%% meake peak mask
fLabels    = unique(fieldsLabel)';
thresholds = nan(size(map));
for i = fLabels(2:end)   
    % thresholds(fieldsLabel == i) = max(map(fieldsLabel == i),[],'omitnan')/2;
    thresholds(fieldsLabel == i) = prctile(map(fieldsLabel == i),68);
end
% generate peak mask and do a bit of cleaning up
tmpMask               = map > thresholds;
tmpMask               = imclose(tmpMask,strel('square',3)); % merge peaks that are too close to each other
% remove pixel bridges 
tmpMask(isnan(map))   = 1;
tmpMask               = ~bwmorph(~tmpMask,'bridge');        
tmpMask(isnan(map))   = 0;
% now we split fields that are only connected on diagonal of 2 pixels
CC                    = bwconncomp(tmpMask,4);      
peakMask              = labelmatrix(CC);
% 
peakStats                                                 = regionprops(peakMask,map,'WeightedCentroid','Area','PixelIdxList','MajorAxisLength','EquivDiameter');
tooSmallFields                                            = [peakStats.Area] < options.minPeakSz;
peakMask(vertcat(peakStats(tooSmallFields).PixelIdxList)) = 0;
peakStats([peakStats.Area] < options.minPeakSz)           = [];


%% debug plot
if options.debugOn
    [peakY,peakX] = ind2sub(size(map),vertcat(peakStats.PixelIdxList)');
    figure;
    subplot(1,2,1);
    scanpix.plot.plotRateMap(map,gca);
    subplot(1,2,2);
    scanpix.plot.plotRateMap(map,gca,'colmap','hcg');
    hold on
    scatter(gca,peakX,peakY,48,'filled','r');
    hold off
end


% switch type
%     case 'place'
%         %% TO DO %% 
%         % need to import old code from the SCAN era
%     case 'grid'
%         % make sAC and grab grid properties
%         sac = scanpix.analysis.spatialCrosscorr(rMap,rMap);
% %         [~,gridProps] = scanpix.analysis.gridprops_v2(sac,'peakMode','point','corrThr',0,'radius','est');
%         [~,gridProps] = scanpix.analysis.gridprops(sac,'legacyMode',true);
%         if isnan(gridProps.fieldSize(1)); peakCoords = []; return; end
%         % gaussian fit of central peak
%         [X, Y] = meshgrid(1:size(sac,1),1:size(sac,2));
%         % 
%         [~, rho] = cart2pol(X(gridProps.centralPeakMask{1}),Y(gridProps.centralPeakMask{1}));
%         pd = fitdist(rho,'Normal');
%         diameter = ceil(2*sqrt(gridProps.fieldSize(1)/pi)); % use field size as proxy for kernel size
%         % covolve rate map with LoG 
%         kernel = fspecial('log',[diameter diameter],pd.sigma);
%         unVisBins = isnan(rMap);
%         rMap(unVisBins) = 0;
%         filtMap = imfilter(rMap,kernel);
%         % reset to NaN
%         rMap(unVisBins) = NaN;
%         filtMap(unVisBins) = 0;
%         % local minima correspond to peak positions
%         bw = imregionalmin(filtMap,8);
%         bw(unVisBins) = 0;
% %         bwL = bwlabel(bw);
%         stats1 = regionprops(bw,rMap,'MaxIntensity','Area','PixelList');
%         % filter fields
%         % rate
%         rateThresh = prctile(rMap(:),75);
%         stats1 = stats1([stats1(:).MaxIntensity] > rateThresh);
%         % overlap - if peaks overlap, only keep the one with higher rate       
%         while ~all(pdist(vertcat(stats1.PixelList)) >= diameter)
%             dist = triu(squareform(pdist(vertcat(stats1.PixelList))));
%             dist(dist==0) = NaN;
%             [r,c] = find(dist<diameter);
%             if stats1(r(1)).MaxIntensity > stats(c(1)).MaxIntensity
%                 stats1(c(1)) = [];
%             else
%                 stats1(r(1)) = [];
%             end
%         end
%         peakCoords = vertcat(stats.PixelList);
% 
%         if prms.debugOn
%             figure;
%             subplot(1,2,1);
%             scanpix.plot.plotRateMap(rMap,gca);
%             axis square;
%             subplot(1,2,2);
%             scanpix.plot.plotRateMap(rMap,gca,'colmap','hcg');
%             axis square;
%             hold on
%             scatter(gca,peakCoords(:,1),peakCoords(:,2),48,'filled','r');
%             hold off
%         end
% end


end

