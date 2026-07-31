function [peakStats, peakMask] = fieldDetect(map,options)
% find peaks in an arbitrary rate map, spatial AC etc. The approach might
% seem very involved but it makes sure that we can find local peaks in e.g.
% a spatial AC that do not have a strong prominence and that peaks that
% seem to bleed into each other are still sorted into different peaks
% package: scanpix.analysis
%
% Usage:
%       [peakStats, peakMask] = scanpix.analysis.gridprops( map );
%       [peakStats, peakMask] = scanpix.analysis.gridprops( map , 'paramName', 'paramValue', .. );
%
% logic:
% - (optional) bin and/or threshold map
% - segment using watershed
% - merge fields (too small) to counter oversegmentation
% - find peaks in basins using a local threshold
% - clean up by merging too close peaks, removing too small peaks
%
% A few notes:
% - generally tends to work best on binned maps for spatial ACs ('binMap' = true, 'thr' = 0, 'thrMode' = 'abs')
% - tends to work best on coarsely binned rate maps ('binMap' = true, 'nBinSteps' = 5)



%%
arguments
    map {mustBeNumeric} 
    options.binMap (1,1) {mustBeNumericOrLogical} = false;
    options.nBinSteps (1,1) {mustBeNumeric} = 11; 
    options.binEdges (1,2) {mustBeNumeric} = [0 max(map(:),[],'omitnan')];
    options.thrMode (1,:) {mustBeMember(options.thrMode,{'abs','rel','none'})} = 'none';
    options.thr {mustBeScalarOrEmpty} = []; 
    options.minWSFieldSz (1,1) {mustBeNumeric} = 16; 
    options.minPeakSz (1,1) {mustBeNumeric} = 8; 
    options.debugOn (1,1) {mustBeNumericOrLogical} = false;
end


%%
if options.binMap
    tmpMap = scanpix.maps.binAnyRMap(map,'nsteps',options.nBinSteps,'cmapEdge',options.binEdges);
else
    tmpMap = map;
end

% mask NaNs for watershed and, if desired, background pixels
switch options.thrMode
    case 'abs'
        tmpMap(map < options.thr | isnan(map)) = -Inf;
        thr                                    = options.thr;
    case 'rel'
        thr                                    =  max(map(:),[],'omitnan') * options.thr;
        tmpMap(map < thr | isnan(map))         = -Inf;
    case 'none'
        tmpMap(isnan(map)) = -Inf;
        thr                = max(map(:),[],'omitnan');

end
% watershed
fieldsLabel = watershed(-tmpMap);
% special case if only a single field is present and the rest of the map is -Inf after thresholding - watershed returns all 1's in that cxase
if all(fieldsLabel(:))
    fieldsLabel = tmpMap ~= -Inf;
end

%% merge fields that are too small into the larger neighbour with the most shared ridge pixels
allLabels   = unique(fieldsLabel(fieldsLabel ~= 0));
fieldSizes  = accumarray(fieldsLabel(fieldsLabel ~= 0),1);
bigLabels   = allLabels(fieldSizes(allLabels) >= options.minWSFieldSz);
smallLabels = allLabels(fieldSizes(allLabels) <  options.minWSFieldSz);

if ~isempty(smallLabels) && ~isempty(bigLabels)
    % process smallest fields first so chained merges (a too-small field
    % bordering another too-small field) resolve into an already-merged
    % big neighbour rather than needing a separate pass
    [~,ord]     = sort(fieldSizes(smallLabels));
    smallLabels = smallLabels(ord);
    centroids   = []; % lazily filled in only if the centroid fallback below is ever needed

    for i = 1:numel(smallLabels)
        lbl           = smallLabels(i);
        fieldMask     = fieldsLabel == lbl;
        % everything within reach of this field's border (ridge pixels and,
        % if the ridge is thin, the neighbouring field's own pixels too).
        % Candidates may include other still-too-small fields - those get
        % their own turn later in the (ascending-size) loop, so the chain
        % still ends up folded into a big field by the time we're done
        dilFieldMask  = quickDilate(quickDilate(fieldMask)); % radius 2 - bridges a multi-pixel-wide ridge
        neighbourLbls = fieldsLabel(dilFieldMask & ~fieldMask & fieldsLabel > 0);

        if ~isempty(neighbourLbls)
            % target = neighbour with the most shared border/ridge pixels
            counts     = accumarray(neighbourLbls,1,[max(allLabels) 1]);
            [~,target] = max(counts);
        else
            % not directly touching anything - fall back to nearest centroid
            if isempty(centroids)
                centroidStats = regionprops(fieldsLabel,'Centroid');
                centroids     = vertcat(centroidStats.Centroid);
            end
            d          = vecnorm(centroids(bigLabels,:) - centroids(lbl,:), 2, 2);
            [~,k]      = min(d);
            target     = bigLabels(k);
        end
        % close the ridge seam between this field and its chosen neighbour,
        % but don't bridge into ridge segments touching a third field. Use a
        % tight (immediate-neighbour, radius 1) test here so a merely-nearby
        % third field doesn't wrongly veto closing a genuine two-field seam
        targetMask = fieldsLabel == target;
        otherMask  = fieldsLabel > 0 & fieldsLabel ~= lbl & fieldsLabel ~= target;
        mergeRidge = quickDilate(fieldMask) & quickDilate(targetMask) & fieldsLabel == 0 & ~quickDilate(otherMask);
        %
        fieldsLabel(fieldMask | mergeRidge) = target;
    end
end

%% meake peak mask
fLabels    = unique(fieldsLabel)';
thresholds = nan(size(map));
for i = fLabels(2:end)   
    thresholds(fieldsLabel == i) = max(thr,prctile(map(fieldsLabel == i),75)); 
end
% generate peak mask and do a bit of cleaning up
tmpMask               = map > thresholds;
tmpMask               = bwareaopen(tmpMask, options.minPeakSz); % remove peaks that are too small
tmpMask               = imclose(tmpMask,strel('square',3));
% tmpMask               = imclose(tmpMask,strel('square',3)); % merge peaks that are too close to each other
% remove pixel bridges 
% tmpMask(isnan(map))   = 1;
% tmpMask               = ~bwmorph(~tmpMask,'bridge');        
tmpMask(isnan(map))   = 0;
% now we split fields that are only connected on diagonal of 2 pixels
CC                    = bwconncomp(tmpMask,4);      
peakMask              = labelmatrix(CC);
% final stats of found peaks 
peakStats                                         = regionprops(peakMask,map,'WeightedCentroid','Area','PixelIdxList','MajorAxisLength','EquivDiameter','Centroid','PixelList','PixelValues');
% remove fields that are too small
tooSmallFields                                    = [peakStats.Area] < options.minPeakSz;
oneBinWideFields                                  = cellfun(@(x) any(x~=0),cellfun(@(x) all(diff(x,[],1)==0),{peakStats.PixelList},'UniformOutput',false));
remInd                                            = tooSmallFields | oneBinWideFields;
peakMask(vertcat(peakStats(remInd).PixelIdxList)) = 0;
peakMask                                          = peakMask > 0;
peakStats(remInd)                                 = [];

% get location of absolute field peak as well
for i = 1 : length(peakStats)
    % Find index of max value
    [~, maxIdx]          = max(peakStats(i).PixelValues,[],'omitnan');
    peakStats(i).peakLoc = peakStats(i).PixelList(maxIdx,:); 
end

%% debug plot
if options.debugOn
    [peakY,peakX] = ind2sub(size(map),vertcat(peakStats.PixelIdxList)');
    figure;
    subplot(1,2,1);
    if options.binMap
        imagesc(gca,map);
        axis square
    else
        scanpix.plot.plotRateMap(map,gca);
    end
    subplot(1,2,2);
    scanpix.plot.plotRateMap(map,gca,'colmap','hcg');
    hold on
    scatter(gca,peakX,peakY,48,'filled','r');
    hold off
end

end

function d = quickDilate(mask)
% 8-connected 1-pixel dilation. Equivalent to imdilate(mask,strel('square',3))
% but avoids imdilate's generic dispatch overhead, which dominates runtime.
p = false(size(mask,1)+2, size(mask,2)+2);
p(2:end-1,2:end-1) = mask;
d = p(1:end-2,1:end-2) | p(1:end-2,2:end-1) | p(1:end-2,3:end) | ...
    p(2:end-1,1:end-2) | p(2:end-1,2:end-1) | p(2:end-1,3:end) | ...
    p(3:end,1:end-2)   | p(3:end,2:end-1)   | p(3:end,3:end);
end

% %% merge fields that are too small into the closest larger one
% fieldMask    = fieldsLabel ~= 0;
% tmpStats     = regionprops(fieldMask,fieldsLabel,'centroid','PixelList','MaxIntensity');
% % sort by size
% [sz,sortInd] = sort( arrayfun(@(x) size(x.PixelList,1), tmpStats) );
% tmpStats     = tmpStats(sortInd);
% % index of fields with too small size
% tooSmallInd  = find(sz < options.minWSFieldSz)';
% 
% % all field boundary pixels
% fieldBorderPix = arrayfun(@(x) bwdist(fieldsLabel==x.MaxIntensity) < 2 & ~(fieldsLabel==x.MaxIntensity),tmpStats, 'UniformOutput',0);
% % fieldBorderPix = arrayfun(@(x) imdilate(fieldsLabel==x.MaxIntensity, strel('square',3)) & ~(fieldsLabel==x.MaxIntensity), tmpStats, 'UniformOutput', 0);
% 
% 
% while ~isempty(tooSmallInd)
%     structInd         = 1:length(tmpStats);
%     % all field centroid distances
%     % centroids         = reshape([tmpStats.Centroid],2,[])';
%     % dists             = squareform(pdist(centroids));
%     % dists(dists == 0) = NaN;
%     % % closest field
%     % [~,minInd]        = min(dists(tooSmallInd(1),:),[],'omitnan');
% 
%     centroids   = reshape([tmpStats.Centroid],2,[])';
%     d           = pdist2(centroids(tooSmallInd(1),:), centroids);
%     d(tooSmallInd(1)) = NaN;
%     [~,minInd]  = min(d,[],'omitnan');
% 
% 
%     % 
%     otherFieldsInd    = ~ismember(structInd,[tooSmallInd(1);minInd]);
%     if ~any(otherFieldsInd)
%         otherFieldsBorders = false(size(tmpMap));
%     else
%         otherFieldsBorders = any(cat(3,fieldBorderPix{otherFieldsInd}),3);
%     end
%     mergeBorderPix    = fieldBorderPix{tooSmallInd(1)} & fieldBorderPix{minInd} & ~otherFieldsBorders;
%     % update fields label
%     fieldsLabel(mergeBorderPix)                                       = tmpStats(minInd).MaxIntensity;
%     fieldsLabel(fieldsLabel == tmpStats(tooSmallInd(1)).MaxIntensity) = tmpStats(minInd).MaxIntensity;
%     %
%     % update structure after merging
%     tmpStats(minInd).PixelList     = [tmpStats(minInd).PixelList; tmpStats(tooSmallInd(1)).PixelList];
%     tmpStats(minInd).Centroid      = mean([tmpStats(minInd).Centroid; tmpStats(tooSmallInd(1)).Centroid],1);
%     tmpStats(tooSmallInd(1))       = [];
%     % update border mask after merging
%     fieldBorderPix{minInd}         = (fieldBorderPix{tooSmallInd(1)} | fieldBorderPix{minInd}) & ~mergeBorderPix; 
%     fieldBorderPix(tooSmallInd(1)) = [];
%     % check for more fields < size thresh
%     tooSmallInd                    = find(arrayfun(@(x) size(x.PixelList,1), tmpStats) < options.minWSFieldSz)';
% end

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

