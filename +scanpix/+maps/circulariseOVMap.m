function ovMapsCirc = circulariseOVMap(obj, trialInd)
% circulariseOVMap - Resample object vector maps (direction x distance) onto
% a cartesian grid centred on the object, i.e. make a 'circular' version
% of the map (as in makeEgoCentricRateMap_v2)
% package: scanpix.maps
%
% Bin centres of the OV map are converted to cartesian coordinates and the
% map is linearly interpolated onto a square grid with a resolution of 1
% distance bin (obj.mapParams.objVect.binSz_dist). Grid points closest to
% an unvisited bin in the OV map, or beyond the map's distance range, are
% set to NaN. The orientation matches the position data / rate maps (x
% along columns, y along rows, i.e. use imagesc to plot), with the object
% at the centre of the grid.
%
% Syntax:
%       ovMapsCirc = scanpix.maps.circulariseOVMap(obj, trialInd)
%
% Inputs:
%    obj         - ephys class object (OV maps for trial need to exist - see obj.addMaps('objVect'))
%    trialInd    - numeric index of trial
%
% Outputs:
%    ovMapsCirc  - nCell x 1 cell array of circular OV maps
%
% see also: scanpix.maps.makeOVMap; scanpix.maps.addMaps
%
% LM 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
end

%% maps & params
if length(obj.maps.OV) < trialInd || isempty(obj.maps.OV{trialInd})
    error('scaNpix::maps::circulariseOVMap:No object vector maps for trial %i. Make some with obj.addMaps(''objVect'') first.',trialInd);
end
ovMaps = obj.maps.OV{trialInd};

% use defaults for any field missing in obj (e.g. params loaded from older files)
prms = scanpix.maps.defaultParamsRateMaps;
prms = prms.objVect;
if isfield(obj.mapParams,'objVect')
    f = fieldnames(obj.mapParams.objVect);
    for i = 1:length(f);   prms.(f{i}) = obj.mapParams.objVect.(f{i});   end
end

mapInd = find(~cellfun('isempty',ovMaps),1);
if isempty(mapInd); ovMapsCirc = cell(size(ovMaps)); return; end
[nDir, nDist] = size(ovMaps{mapInd});
if nDir ~= round(360 / prms.binSz_dir)
    error('scaNpix::maps::circulariseOVMap:Map size (%i direction bins) doesn''t match binSz_dir (%g deg). Did you change the params after making the maps?',nDir,prms.binSz_dir);
end

%% polar bin centres -> cartesian (in units of distance bins, object at 0,0)
thetaC          = ((1:nDir)' - 0.5) .* prms.binSz_dir .* pi/180;
rMin            = prms.minDist / prms.binSz_dist;  % inner edge of map
rC              = rMin + (1:nDist) - 0.5;
rMax            = rMin + nDist;                    % outer edge of map
[R, TH]         = meshgrid(rC, thetaC);            % nDir x nDist, same as map
[binX, binY]    = pol2cart(TH, R);

% resampling grid
nR              = ceil(rMax);
[gridX, gridY]  = meshgrid(-nR:nR);
gridR           = sqrt(gridX.^2 + gridY.^2);
outOfRange      = gridR > rMax | gridR < rMin;

%% interpolate
ovMapsCirc = cell(size(ovMaps));
F          = [];
for i = 1:length(ovMaps)
    if isempty(ovMaps{i}); continue; end
    valid = ~isnan(ovMaps{i});
    if nnz(valid) < 3
        ovMapsCirc{i} = nan(size(gridX));
        continue
    end
    % unvisited bins are identical for all cells in a trial, so only (re)build interpolants if necessary
    if isempty(F) || ~isequal(valid, validF)
        validF  = valid;
        F       = scatteredInterpolant(binX(valid), binY(valid), ovMaps{i}(valid), 'linear', 'none');
        % keep unvisited parts of the map unvisited (interpolation would otherwise fill them)
        Fvis    = scatteredInterpolant(binX(:), binY(:), double(valid(:)), 'nearest', 'none');
        unVis   = reshape(Fvis(gridX(:), gridY(:)), size(gridX)) ~= 1 | outOfRange;
    else
        F.Values = ovMaps{i}(valid);
    end
    ovMapsCirc{i}        = reshape(F(gridX(:), gridY(:)), size(gridX));
    ovMapsCirc{i}(unVis) = NaN;
end

end
