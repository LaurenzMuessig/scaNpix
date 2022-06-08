function objMap = makeOVMap(spikeTimes,xy,sampleT,objPos,ppm,varargin)


%% params

%% params
prms.binSz_dist = 2.5; % in cm;  2cm in Høydal et al (2019)
prms.minDist    = 5;   % in cm;
prms.binSz_dir  = 10;  % in degrees;  5deg in Høydal et al (2019)
prms.posFs      = 50;  % sample rate
% smoothing
prms.smKernelSz_OV = 5;
prms.smSigma_OV    = 2;
prms.showWaitBar   = false;

% prms.debugOn    = 0;

%% parse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(spikeTimes)
    spikeTimes = {spikeTimes};
end

if length(prms.smKernelSz_OV) == 1
    prms.smKernelSz_OV = [prms.smKernelSz_OV prms.smKernelSz_OV];
end

% this is a temp hack to deal with multiple obj trials - for now use first
% in list as hard-coded
if length(objPos) > 2
    warning('Coordinates for several objects supplied. Will use first in list as reference. Multi object detection is not yet supported!');  
    objPos = objPos(1:2);
end

%% PREPROCESS
% objPos               = objPos - min(xy);
% xy                   = ceil( bsxfun(@minus,xy, min(xy)) + eps ); % make min pos data = [1,1]
binSzDir             = prms.binSz_dir * pi/180;
% binSzDist            = ceil(prms.binSz_dist * (ppm/100)); % in camera pixel

%% MAKE OCCUPANCY MAP
% distances to object
dist                 = sqrt( (xy(:,1) - objPos(1)).^2 + (xy(:,2) - objPos(2)).^2 ) ./ (ppm/100); % all distances to obj in cm
notValidInd          = dist < prms.minDist;  
xy(notValidInd,:)    = NaN;
dist                 = ceil( dist ./ prms.binSz_dist ); % binned
minDistBInned        = ceil( prms.minDist ./ prms.binSz_dist );
% angles to object      
theta                = mod(atan2(xy(:,2)-objPos(2),xy(:,1)-objPos(1)), 2*pi); % all angles to obj in degrees %% IS THIS RIGHT??

% theta(theta < 0)     = theta(theta < 0) + 2*pi; % 0:360
theta                = ceil( theta ./ binSzDir ); % binned
% occupancy map 
occMap               = accumarray([theta(~isnan(xy(:,1))) dist(~isnan(xy(:,1))) ],1,[nanmax(theta(:)) nanmax(dist(:))]) ./ prms.posFs; 
occMap               = occMap(:,minDistBInned:end);
% loop over cells to make rate maps
[ spkMaps, objMap ]   = deal(cell(length(spikeTimes),1));

if prms.showWaitBar; hWait = waitbar(0); end

for i = 1:length(spikeTimes)
    if isempty(spikeTimes{i})
        [objMap{i}, spkMaps{i}] = deal(zeros(size(occMap)));
        continue
    end
    if isempty(sampleT)
        spkPosBinInd = ceil(spikeTimes{i} .* prms.posFs ); 
    else
%         [~, spkPosBinInd] = arrayfun(@(x) min(abs(sampleT - x)), spikeTimes{i}, 'UniformOutput', 0);
%         spkPosBinInd = cell2mat(spkPosBinInd);
        [~, spkPosBinInd] = min(abs(bsxfun(@minus, sampleTimes, spikeTimes{i}.')), [], 1);
    end
    spkBinnedDist    = dist(spkPosBinInd); 
    spkBinnedTheta   = theta(spkPosBinInd); 
    nanInd           = isnan(spkBinnedTheta) | isnan(spkBinnedDist);
    spkMaps{i}       = accumarray([spkBinnedTheta(~nanInd) spkBinnedDist(~nanInd)],1,[nanmax(theta(:)) nanmax(dist(:))]);
    spkMaps{i}       = spkMaps{i}(:,minDistBInned:end);
    objMapRaw        = spkMaps{i} ./ occMap;
    objMapRaw(occMap==0) = 0;
    % smooth map
    % as map is circular-linear we need to padd both dims differently
    objMapRawPadded  = padarray(objMapRaw,prms.smKernelSz_OV(1),'circular');
    objMapRawPadded  = padarray(objMapRawPadded,[0 prms.smKernelSz_OV(2)],0);  
    temp             = imgaussfilt(objMapRawPadded,prms.smSigma_OV,'filtersize', prms.smKernelSz_OV);  % smooth map
    objMap{i}        = temp(prms.smKernelSz_OV(1)+1:end-prms.smKernelSz_OV(1),prms.smKernelSz_OV(2)+1:end-prms.smKernelSz_OV(2)); % remove padding
    objMap{i}(occMap==0) = nan;

    if prms.showWaitBar; waitbar(i/length(spikeTimes),hWait,sprintf('Making those Object Vector Maps... %i/%i done.',i,length(spikeTimes))); end
    
%     if prms.debugOn
%         figure; 
%         subplot(1,2,1);
%         rMap = makeRateMaps(spikeTimes, xy, sampleT, ppm);
%         [rMapBinned, cMapBinned] = binAnyRMap(rMap{1},'jet',11,[1 1 1]);
%         imagesc(rMapBinned); colormap(gca,cMapBinned); axis square; axis off
%         hold on
%         objPosBins = objPos ./ 10;
%         scatter(objPosBins(2),objPosBins(1),108,'rx','linewidth',4);
%         hold off
%         subplot(1,2,2);
%         imagesc(objMap{i},[0 nanmax(objMap{i}(:))]); axis square
%         set(gca,'xtick',5:5:size(objMap{i},2),'XTickLabel',(5:5:size(objMap{i},2)).*prms.binSz_dist,'ytick',0.5:5:60.5,'YTickLabel',0:30:360);
%     end
end

if prms.showWaitBar; close(hWait); end

end

