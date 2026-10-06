function [objMap, occMap, spkMaps] = makeOVMap(obj, trialInd, options)
% makeOVMap - Make object vector maps, i.e. firing rate as a function of
% direction and distance from an object (Hoydal et al. 2019, Nature)
% package: scanpix.maps
%
% Direction is allocentric, from object to animal, with 0 = positive x
% direction in camera coordinates. Distance is measured from the object
% centre. Object coordinates (trialMetaData.objectPos) are expected in raw
% camera pixels and are converted here to match the processed position
% data (scanpix.maps.rawToFitCoords). Map params are taken from
% obj.mapParams.objVect.
%
% NOTE: scanpix.maps.scalePosition does not only shift the positions, it
% also rescales them. The full transform is stored in
% trialMetaData.PosIsFitToEnv{3} since 2026-10 - for data loaded before
% that only the shift is known, so the object position will be slightly off
% for rectangular envs (~1-2% stretch, i.e. up to ~1cm) and can be off
% substantially for circular envs (rawToFitCoords warns). Reload to fix.
%
% Syntax:
%       objMap = scanpix.maps.makeOVMap(obj, trialInd)
%       [objMap, occMap, spkMaps] = scanpix.maps.makeOVMap(obj, trialInd, Name-Value comma separated list)
%
% Inputs:
%    obj          - ephys class object
%    trialInd     - numeric index of trial
%    options      - name-value: 'addPosFilter' (logical nPosx1, positions to exclude)
%                               'cellInd' (logical/numeric index of cells to make maps for)
%
% Outputs:
%   objMap       - nCell x 1 cell array of smoothed object vector maps (rows: direction, cols: distance)
%   occMap       - occupancy map (s)
%   spkMaps      - nCell x 1 cell array of raw spike count maps
%
% see also: scanpix.maps.addMaps; scanpix.maps.makeRateMaps; scanpix.maps.rotatePosition
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
    options.addPosFilter {mustBeNumericOrLogical} = false(size(obj.posData.XY{trialInd},1),1);
    options.cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);
end

%% params
% use defaults for any field missing in obj (e.g. params loaded from older files)
prms = scanpix.maps.defaultParamsRateMaps;
prms = prms.objVect;
if isfield(obj.mapParams,'objVect')
    f = fieldnames(obj.mapParams.objVect);
    for i = 1:length(f);   prms.(f{i}) = obj.mapParams.objVect.(f{i});   end
end
prms.posFs = obj.trialMetaData(trialInd).posFs;

% only relevant for npix data - check for interpolated Fs
if strcmp(obj.type,'npix') && isKey(obj.params,'InterpPos2PosFs') && obj.params('InterpPos2PosFs')
    sampleT = [];
else
    sampleT = obj.spikeData.sampleT{trialInd};
end

% data from object
xy                         = obj.posData.XY{trialInd};
xy(options.addPosFilter,:) = NaN;
ppm                        = obj.trialMetaData(trialInd).ppm;
spikeTimes                 = obj.spikeData.spk_Times{trialInd}(options.cellInd);

%% object position
if ~isfield(obj.trialMetaData(trialInd),'objectPos') || numel(obj.trialMetaData(trialInd).objectPos) < 2
    error('scaNpix::maps::makeOVMap:No object position (trialMetaData.objectPos) for trial %i.',trialInd);
end
% stored as X1Y1,X2Y2,... in raw camera pixels
objPos = obj.trialMetaData(trialInd).objectPos;
nObj   = floor(numel(objPos)/2);
objPos = [reshape(objPos(1:2:2*nObj),[],1), reshape(objPos(2:2:2*nObj),[],1)];
% this is a temp hack to deal with multiple obj trials - for now use first in list as hard-coded
if nObj > 1
    warning('scaNpix::maps::makeOVMap:Coordinates for several objects supplied. Will use first in list as reference. Multi object detection is not yet supported!');
end
% convert to the frame of the position data (scaling to common ppm and fit to environment)
objPos = scanpix.maps.rawToFitCoords(obj, trialInd, objPos(1,:));

%% speed filter
if prms.speedFilterFlagOVMaps
    speedFilter       = obj.posData.speed{trialInd} <= prms.speedFilterLimitLow | obj.posData.speed{trialInd} > prms.speedFilterLimitHigh;
    xy(speedFilter,:) = NaN;
end

if all(isnan(xy(:,1)))
    warning('scaNpix::maps::makeOVMap:No valid position samples left for trial %i (check speed filter). No maps generated.',trialInd);
    [objMap, spkMaps] = deal(cell(length(spikeTimes),1));
    occMap            = [];
    return
end

%% BIN DISTANCE & DIRECTION
nDirBins = round(360 / prms.binSz_dir);
if abs(nDirBins * prms.binSz_dir - 360) > 1e-9
    error('scaNpix::maps::makeOVMap:binSz_dir (%g deg) needs to divide 360.',prms.binSz_dir);
end
% distances (cm) and angles (rad, object -> animal, 0 = +x) to object
dist  = sqrt( (xy(:,1) - objPos(1)).^2 + (xy(:,2) - objPos(2)).^2 ) ./ (ppm/100);
theta = mod(atan2(xy(:,2)-objPos(2), xy(:,1)-objPos(1)), 2*pi);
% distance bins start at minDist; maxDist sets a fixed map size (use this when comparing maps across trials)
if isempty(prms.maxDist)
    nDistBins = floor( (max(dist,[],'omitnan') - prms.minDist) / prms.binSz_dist ) + 1;
else
    nDistBins = ceil( (prms.maxDist - prms.minDist) / prms.binSz_dist );
end
% bin index = floor(x/binSz)+1, so x=0 falls in bin 1
distBin  = floor( (dist - prms.minDist) ./ prms.binSz_dist ) + 1;
thetaBin = min( floor( theta ./ (prms.binSz_dir*pi/180) ) + 1, nDirBins); % min() guards against rounding at 2*pi
valid    = ~isnan(distBin) & ~isnan(thetaBin) & distBin >= 1 & distBin <= nDistBins;
distBin(~valid)  = NaN;
thetaBin(~valid) = NaN;
mapSz    = [nDirBins nDistBins];

%% OCCUPANCY MAP
occMap   = accumarray([thetaBin(valid) distBin(valid)], 1, mapSz) ./ prms.posFs;
unVisPos = occMap == 0;

%% SMOOTHING KERNEL
% Gaussian, SD in bins (scalar or [dir dist]); kernel size defaults to 2*ceil(2*SD)+1 (as imgaussfilt)
smSigma = prms.smSigma_OV .* [1 1];
if isempty(prms.smKernelSz_OV)
    kSz = 2*ceil(2*smSigma) + 1;
else
    kSz = prms.smKernelSz_OV .* [1 1];
    kSz = kSz + (mod(kSz,2) == 0); % needs to be odd
end
halfK  = (kSz - 1) / 2;
gDir   = exp( -(-halfK(1):halfK(1)).^2 ./ (2*smSigma(1)^2) )';
gDist  = exp( -(-halfK(2):halfK(2)).^2 ./ (2*smSigma(2)^2) );
kernel = gDir * gDist;
kernel = kernel ./ sum(kernel(:));
% spike and occupancy maps are smoothed separately and then divided (rather than smoothing the rate map as in
% Hoydal et al.), as polar bins close to the object are tiny and their raw rates are very noisy. Zero padding along
% the distance axis is the same for both maps, so cancels in the ratio
occMap_sm = smoothCircLin(occMap, kernel, halfK);

%% SPIKE TIMES -> POS SAMPLES
nPos = size(xy,1);
if ~isempty(sampleT)
    sampleT = sampleT(:);
    % nearest pos sample for each spike - spikes outside the pos sample range return NaN and are discarded
    posInd = @(t) interp1(sampleT, (1:length(sampleT))', t(:), 'nearest');
else
    % pos sample k covers [(k-1)/posFs, k/posFs)
    posInd = @(t) floor(t(:) .* prms.posFs) + 1;
end

%% RATE MAPS
[ spkMaps, objMap ]  = deal(cell(length(spikeTimes),1));

if prms.showWaitBar; hWait = waitbar(0); end

for i = 1:length(spikeTimes)
    spkPosBinInd = posInd(spikeTimes{i});
    spkPosBinInd = spkPosBinInd(spkPosBinInd >= 1 & spkPosBinInd <= nPos);
    spkTheta     = thetaBin(spkPosBinInd);
    spkDist      = distBin(spkPosBinInd);
    ok           = ~isnan(spkTheta);
    spkMaps{i}   = accumarray([spkTheta(ok) spkDist(ok)], 1, mapSz);
    % smoothed spikes / smoothed occupancy (circular in direction)
    objMap{i}           = smoothCircLin(spkMaps{i}, kernel, halfK) ./ occMap_sm;
    objMap{i}(unVisPos) = NaN;

    if prms.showWaitBar; waitbar(i/length(spikeTimes),hWait,sprintf('Making those Object Vector Maps... %i/%i done.',i,length(spikeTimes))); end
end

if prms.showWaitBar; close(hWait); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mapSm = smoothCircLin(map, kernel, halfK)
% smooth a direction x distance map: circular padding along direction (rows), zero padding along distance (cols)
mapPadded = [map(end-halfK(1)+1:end,:); map; map(1:halfK(1),:)];
mapPadded = [zeros(size(mapPadded,1),halfK(2)), mapPadded, zeros(size(mapPadded,1),halfK(2))];
mapSm     = conv2(mapPadded, kernel, 'valid');
end
