function [gridness, Props] = gridprops_v2(autoCorr, options)
% gridprops_v2 - grid properties from the spatial autocorrelogram (SAC), with the
% inner six peaks found via the troughs around the central peak (gridPeaksFromTroughs).
%
%       [gridness, Props] = scanpix.analysis.gridprops_v2(autoCorr);
%       [gridness, Props] = scanpix.analysis.gridprops_v2(autoCorr, 'minContrast', 0.2);
%       prms   = scanpix.maps.defaultParamsRateMaps;
%       gpArgs = namedargs2cell(prms.gridPropsV2);
%       [gridness, Props] = scanpix.analysis.gridprops_v2(autoCorr, gpArgs{:});
%
% package: scanpix.analysis   (plot with scanpix.plot.plotGridProps_v2)
%
% Peaks: regional maxima whose ray trough sits at |P|/2, combined into the lattice
% (a1, a2, a1+-a2) with the best peak-vs-hole contrast. If no lattice passes the
% rejection rules (min peak prominence, min lattice contrast) no peaks are returned
% (wavelength/orientation/ellipse = NaN), but gridness is still computed.
%
% Gridness (Moser style, as gridprops): correlate the SAC with its 30..150 deg rotations
% within an annulus; inner edge = trough ring (edge of the central peak), outer edge
% grown iteratively from just beyond the peak ring to the SAC edge; gridness = max over
% outer radii of min(r60,r120) - max(r30,r90,r150). The annulus never starts narrower
% than 'minAnnWidth' x inner radius, so it always contains the peak ring.
%
% Elliptical gridness: SAC is compressed along the major axis of the ellipse through the
% six peaks (centre fixed at origin; censored peaks replaced by 2x trough vector), so the
% peaks lie on a circle, then gridness is computed as above. The correction is used if the
% ellipse ratio is <= maxEllRatio (1.5), or <= maxEllRatioSupported (2.5) when all three axes
% have a detected peak that is not cut off by the SAC edge. Otherwise gridnessEll = gridness
% (ellCorrected = false); elongated fits on weakly supported lattices are mostly band-like or
% irregular SACs. The stretched SAC (acReg) and its score (gridnessEllStretched) are always
% returned. No lattice -> NaN.
%
% Units: all distances in SAC bins; coordinates x,y in image coords (y down).
%
% Props fields (same names as gridprops where possible):
%   gridness, gridnessEll, wavelength (per axis: peak distance, or 2*trough if the peak
%   is cut off by the SAC edge; mean of 3), wavelengthFull [3x1], wavelengthPeak,
%   wavelengthTrough (2 x mean trough), orientation, orientationFull, offset, offsetFull,
%   peakCoords [7x2] (row 1 = centre), closestPeaksCoord [6x2], censored [3x1],
%   troughCoords [6x2], ellAbScale [a b] (semi axes), ellOrient (rad), ellRatio,
%   acReg (regularised SAC), centralPeakMask, gridMask, gridMaskEll, annulus [rIn rOut],
%   annulusEll, latticeFound, contrast, method, ellCorrected, peakCoordsEll / troughCoordsEll
%   (peaks/troughs in the ellipse corrected SAC), wavelengthEll, wavelengthFullEll, orientationEll,
%   orientationFullEll, offsetEll, offsetFullEll, closestPeaksCoordEll (properties of the corrected
%   SAC), gridnessEllStretched (gridness of the stretched SAC, also when ratio > maxEllRatio)
%
% LM / Claude 2026

arguments
    autoCorr {mustBeNumeric}
    options.cRange (1,2) {mustBeNumeric} = [0.75 1.2];  % allowed |peak| / (2*trough) along a ray (edge peaks read short)
    options.band (1,2) {mustBeNumeric} = [1.2 3.2];     % candidate peaks searched in this range x central trough radius
    options.minSep (1,1) {mustBeNumeric} = 35;          % min angle (deg) between two lattice axes
    options.wedge (1,1) {mustBeNumeric} = 3;            % angular average (deg) for ray profiles
    options.edgeTol (1,1) {mustBeNumeric} = 1;          % peaks within this many bins of the NaN edge are 'censored'
    options.require3 (1,1) {mustBeNumericOrLogical} = false; % reject lattice if third axis has no detected peak
    options.minContrast (1,1) {mustBeNumeric} = 0.25;  % lattice rejection threshold (peak - hole SAC); 96% of grid cells (G>0.5) pass, ~45% of non-grid cells
    options.minProm (1,1) {mustBeNumeric} = 0.1;       % peak must exceed its ray trough by this
    options.minAnnWidth (1,1) {mustBeNumeric} = 1.5;   % first outer radius >= rIn * (1 + minAnnWidth)
    options.maxEllRatio (1,1) {mustBeNumeric} = 1.5;   % ellipse correction always used up to this major/minor ratio
    options.maxEllRatioSupported (1,1) {mustBeNumeric} = 2.5; % ... and up to this ratio if all 3 axes have a detected, uncut peak
    options.fastGridness (1,1) {mustBeNumericOrLogical} = true; % precomputed rotation matrices + one-pass annulus correlations (identical scores, several x faster)
    options.commonMask (1,:) char {mustBeMember(options.commonMask, {'none', 'ell', 'all'})} = 'none'; % compare all rotations on the same bins (valid in every rotation): never / ellipse-corrected SAC only / both
    options.ellMode (1,:) char {mustBeMember(options.ellMode, {'compress', 'stretch', 'area', 'v1'})} = 'compress'; % how the ellipse is made circular ('v1' = gridprops: stretch minor axis, missing data -> 0)
end

gridness = NaN;
sz = size(autoCorr); c = ceil(sz/2);
Props = struct('gridness', NaN, 'gridnessEll', NaN, 'wavelength', NaN, 'wavelengthFull', nan(3,1), ...
    'wavelengthPeak', NaN, 'wavelengthTrough', NaN, 'orientation', NaN, 'orientationFull', nan(3,1), ...
    'offset', NaN, 'offsetFull', nan(3,1), 'peakCoords', nan(7,2), 'closestPeaksCoord', nan(6,2), ...
    'censored', nan(3,1), 'troughCoords', nan(6,2), 'ellAbScale', nan(1,2), 'ellOrient', NaN, 'ellRatio', NaN, ...
    'acReg', nan(sz), 'centralPeakMask', false(sz), 'gridMask', false(sz), 'gridMaskEll', false(sz), ...
    'annulus', nan(1,2), 'annulusEll', nan(1,2), 'latticeFound', false, 'contrast', NaN, 'method', '', 'ellCorrected', false, ...
    'peakCoordsEll', nan(6,2), 'troughCoordsEll', nan(6,2), 'gridnessEllStretched', NaN, ...
    'wavelengthEll', NaN, 'wavelengthFullEll', nan(3,1), 'orientationEll', NaN, 'orientationFullEll', nan(3,1), ...
    'offsetEll', NaN, 'offsetFullEll', nan(3,1), 'closestPeaksCoordEll', nan(6,2), 'found3', false, 'nCensored', NaN);
if isempty(autoCorr) || all(isnan(autoCorr(:))); return; end

% --- peaks via troughs -------------------------------------------------------
G = scanpix.analysis.gridPeaksFromTroughs(autoCorr, 'fallback', false, 'minContrast', options.minContrast, 'minProm', options.minProm, ...
    'cRange', options.cRange, 'band', options.band, 'minSep', options.minSep, 'wedge', options.wedge, 'edgeTol', options.edgeTol, ...
    'require3', options.require3);
Props.contrast = G.contrast; Props.method = G.method; Props.found3 = G.found3; Props.nCensored = G.nCensored;
found = strcmp(G.method, 'lattice');
Props.latticeFound = found;
minSide = min(sz);

if found
    pXY  = G.peakXY - [c(2) c(1)];               % 6x2 relative, rows 1-3 = axes, 4-6 = mirror
    tXY  = G.troughXY - [c(2) c(1)];
    cens = G.censored(1:3);
    tR   = G.troughR(1:3);
    pR   = G.peakR(1:3);
    Props.wavelengthPeak   = mean(pR(~cens));
    Props.wavelengthTrough = 2*mean(tR, 'omitnan');
    Props.censored         = cens;
    Props.peakCoords       = [c(2) c(1); G.peakXY];
    Props.troughCoords     = G.troughXY;
    Gm = latticeGeometry(pXY, tXY, cens, c);
    for f = fieldnames(Gm)'; Props.(f{1}) = Gm.(f{1}); end
    % annulus geometry
    rIn  = mean(sqrt(sum(tXY.^2, 2)), 'omitnan');
    vPk  = pXY(1:3,:); vPk(cens,:) = 2*tXY(cens,:);
    rPk  = max(sqrt(sum(vPk.^2, 2)));
else
    % no lattice: central peak from the angle averaged profile, annulus out to the edge
    rIn = G.r0;
    if isnan(rIn) || rIn > minSide/4; rIn = minSide/4; end
    rPk = 2*rIn;
end
if options.fastGridness; agFun = @annulusGridnessFast; else; agFun = @annulusGridness; end
[gridness, mask, rOut, cMask] = agFun(autoCorr, rIn, rPk, options.minAnnWidth, strcmp(options.commonMask, 'all'));
Props.gridness = gridness; Props.gridMask = mask; Props.centralPeakMask = cMask; Props.annulus = [rIn rOut];

% --- elliptical gridness -------------------------------------------------------
if found && isfinite(G.ellPeak.a)
    E = G.ellPeak;
    Props.ellAbScale = [E.a E.b]; Props.ellOrient = E.orient; Props.ellRatio = E.ratio;
    Rm = [cos(E.orient) -sin(E.orient); sin(E.orient) cos(E.orient)];
    % map the ellipse onto a circle: 'compress' major axis (peaks -> radius b; samples beyond the SAC
    % edge -> NaN), 'stretch' minor axis (peaks -> radius a; no new NaNs but ring moves outward),
    % 'area' both by sqrt(a/b) (peaks -> radius sqrt(ab))
    switch options.ellMode
        case 'compress'; sc = [E.b/E.a 1];
        case {'stretch', 'v1'}; sc = [1 E.a/E.b];
        case 'area';     sc = [sqrt(E.b/E.a) sqrt(E.a/E.b)];
    end
    M  = Rm * diag(sc) * Rm';
    Mi = Rm * diag(1./sc) * Rm';
    [X, Y] = meshgrid((1:sz(2)) - c(2), (1:sz(1)) - c(1));
    src = Mi * [X(:)'; Y(:)'];
    if strcmp(options.ellMode, 'v1')
        % as gridprops' regularise_eliptic_grid: NaN -> 0 before resampling, outside -> 0, NaNs not restored
        ac0 = autoCorr; ac0(isnan(ac0)) = 0;
        acReg = reshape(interp2(ac0, src(1,:) + c(2), src(2,:) + c(1), 'linear', 0), sz);
    else
        acReg = reshape(interp2(autoCorr, src(1,:) + c(2), src(2,:) + c(1), 'linear'), sz);
    end
    Props.acReg = acReg;
    pX2 = (M * pXY')'; tX2 = (M * tXY')'; vP2 = (M * vPk')';
    Props.peakCoordsEll   = pX2 + [c(2) c(1)];
    Props.troughCoordsEll = tX2 + [c(2) c(1)];
    Ge = latticeGeometry(pX2, tX2, cens, c);
    for f = fieldnames(Ge)'; Props.([f{1} 'Ell']) = Ge.(f{1}); end
    rIn2 = mean(sqrt(sum(tX2.^2, 2)), 'omitnan');
    rPk2 = max(sqrt(sum(vP2.^2, 2)));
    [ge, maskE, rOutE] = agFun(acReg, rIn2, rPk2, options.minAnnWidth, ~strcmp(options.commonMask, 'none'));
    Props.gridnessEllStretched = ge; Props.gridMaskEll = maskE; Props.annulusEll = [rIn2 rOutE];
    wellSupported = G.found3 && G.nCensored == 0;
    if E.ratio <= options.maxEllRatio || (E.ratio <= options.maxEllRatioSupported && wellSupported)
        Props.gridnessEll = ge; Props.ellCorrected = true;
    else
        % elongated fits on a poorly supported lattice (inferred third axis / peaks cut off by the
        % SAC edge) mostly come from band-like or irregular SACs: stretching them manufactures
        % hexagonal structure, so the uncorrected gridness is reported
        % (the stretched SAC and its gridness are still returned: acReg, gridnessEllStretched)
        Props.gridnessEll = gridness;
    end
end
end

% =============================================================================
function Gm = latticeGeometry(pXY, tXY, cens, c)
% scale / orientation / offset from the 6 peaks (relative coords, rows 1-3 = axes, 4-6 mirror).
% Scale per axis = peak distance, or 2 x trough distance if the peak is cut off by the SAC edge.
pR = sqrt(sum(pXY(1:3,:).^2, 2)); tR = sqrt(sum(tXY(1:3,:).^2, 2));
wl = pR; wl(cens) = 2*tR(cens);
Gm.wavelengthFull = wl;
Gm.wavelength     = mean(wl, 'omitnan');
% orientation / offset as gridprops (3 peaks with smallest |angle|, y axis flipped)
o = atan2(pXY(:,2), pXY(:,1));
[~, s] = sort(abs(o));
ori = -o(s(1:3));
Gm.orientationFull   = ori;
Gm.orientation       = angle(mean(exp(1i*ori)));
Gm.offsetFull        = pi/4 - abs(mod(ori, pi/2) - pi/4);
Gm.offset            = min(Gm.offsetFull);
Gm.closestPeaksCoord = pXY(s,:) + [c(2) c(1)];      % sub-bin (gridprops v1 rounds with ceil)
end

% =============================================================================
function [gridness, mask, rOutBest, cMask] = annulusGridnessFast(sac, rIn, rPk, minW, common)
% Same result as annulusGridness, faster:
%  - rotations (bilinear, 'crop', as imrotate) are sparse matrices built once per SAC size
%    and cached (identical to imrotate to ~1e-13)
%  - correlations for all outer radii in one pass: bins sorted by distance from the centre,
%    running sums of a, b, a^2, b^2, ab -> Pearson r for every radius
if nargin < 5; common = false; end
sz = size(sac); c = ceil(sz/2); minSide = min(sz);
R = rotationCache(sz);
cMask = R.distMap < rIn;
rEdge = ceil(minSide/2);
first = max(rPk + rIn/2, rIn*(1 + minW));
radii = ceil(min(first, rEdge)):rEdge;
nA    = numel(R.angles);
rot   = nan(numel(sac), nA);
s0 = sac(:); isn = isnan(s0); s0(isn) = 0;
for j = 1:nA
    tmp = R.W{j} * s0;
    if any(isn); tmp(R.S{j} * double(isn) > 0) = NaN; end   % NaN if any weighted neighbour is NaN
    tmp(~R.inMask{j}) = NaN;
    rot(:,j) = tmp;
end
valid = ~isnan(sac(:)) & ~cMask(:);
if common; valid = valid & all(~isnan(rot), 2); end
annCorr = nan(numel(radii), nA);
d = R.distMap(:);
for j = 1:nA
    m  = valid & ~isnan(rot(:,j));
    [ds, o] = sort(d(m));
    a  = sac(m); a = a(o);  b = rot(m, j); b = b(o);
    Sa = cumsum(a); Sb = cumsum(b); Saa = cumsum(a.^2); Sbb = cumsum(b.^2); Sab = cumsum(a.*b);
    for i = 1:numel(radii)
        n = sum(ds < radii(i));
        if n < 10; continue; end
        num = n*Sab(n) - Sa(n)*Sb(n);
        den = sqrt((n*Saa(n) - Sa(n)^2) * (n*Sbb(n) - Sb(n)^2));
        annCorr(i,j) = num / den;
    end
end
[gridness, iBest] = max(min(annCorr(:,1:2), [], 2) - max(annCorr(:,3:end), [], 2));
if isempty(gridness) || isnan(gridness); gridness = NaN; mask = false(sz); rOutBest = NaN; return; end
rOutBest = radii(iBest);
mask = ~cMask & R.distMap < rOutBest;
end

% =============================================================================
function R = rotationCache(sz)
% sparse bilinear rotation matrices for the gridness angles, per SAC size (persistent cache)
persistent cache
if isempty(cache); cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
key = sprintf('%dx%d', sz);
if isKey(cache, key); R = cache(key); return; end
R.angles = [60 120 30 90 150];
cc = (sz + 1)/2; c = ceil(sz/2);
[C, Rr] = meshgrid(1:sz(2), 1:sz(1));
R.distMap = sqrt((Rr - c(1)).^2 + (C - c(2)).^2);
x = C(:) - cc(2); y = Rr(:) - cc(1);
n = prod(sz);
for j = 1:numel(R.angles)
    th = -R.angles(j)*pi/180;
    xs =  cos(th)*x + sin(th)*y + cc(2);           % source col for each output pixel
    ys = -sin(th)*x + cos(th)*y + cc(1);           % source row
    x0 = floor(xs); y0 = floor(ys); fx = xs - x0; fy = ys - y0;
    I = []; J = []; V = [];
    nb = {[0 0], [1 0], [0 1], [1 1]};
    for q = 1:4
        xi = x0 + nb{q}(1); yi = y0 + nb{q}(2);
        w  = (nb{q}(1)*fx + (1-nb{q}(1))*(1-fx)) .* (nb{q}(2)*fy + (1-nb{q}(2))*(1-fy));
        ok = xi >= 1 & xi <= sz(2) & yi >= 1 & yi <= sz(1) & w > 1e-9; % outside = zero padding; negligible weights dropped
        I = [I; find(ok)]; J = [J; sub2ind(sz, yi(ok), xi(ok))]; V = [V; w(ok)]; %#ok<AGROW>
    end
    R.W{j} = sparse(I, J, V, n, n);                % interpolation weights
    R.S{j} = sparse(I, J, 1, n, n);                % support for NaN propagation (output NaN if any weighted neighbour is NaN;
                                                   % imrotate differs only in floating point borderline cases)
    m = imrotate(true(sz), R.angles(j), 'nearest', 'crop');
    R.inMask{j} = m(:);
end
cache(key) = R;
end

% =============================================================================
function [gridness, mask, rOutBest, cMask] = annulusGridness(sac, rIn, rPk, minW, common)
% iterative annulus gridness. Inner edge rIn (trough ring); outer radii from
% max(rPk + rIn/2, rIn*(1+minW)) to the SAC edge. common = true: all rotations are
% compared on the same bins (valid in the SAC and in every rotated copy).
if nargin < 5; common = false; end
sz = size(sac); c = ceil(sz/2); minSide = min(sz);
[cols, rows] = meshgrid(1:sz(2), 1:sz(1));
distMap = sqrt((rows - c(1)).^2 + (cols - c(2)).^2);
cMask   = distMap < rIn;
rotAngle = [60 120 30 90 150];
rot = nan([sz numel(rotAngle)]);
for i = 1:numel(rotAngle)
    tmp = imrotate(sac, rotAngle(i), 'bilinear', 'crop');
    tmp(~imrotate(true(sz), rotAngle(i), 'nearest', 'crop')) = NaN;
    rot(:,:,i) = tmp;
end
rEdge = ceil(minSide/2);
first = max(rPk + rIn/2, rIn*(1 + minW));
radii = ceil(min(first, rEdge)):rEdge;
annCorr = nan(numel(radii), numel(rotAngle));
valid = ~isnan(sac);
if common; valid = valid & all(~isnan(rot), 3); end
for i = 1:numel(radii)
    ann = ~cMask & distMap < radii(i);
    for j = 1:numel(rotAngle)
        m = ann & valid & ~isnan(rot(:,:,j));
        if nnz(m) < 10; continue; end
        r2 = rot(:,:,j);
        annCorr(i,j) = corr2(sac(m), r2(m));
    end
end
[gridness, iBest] = max(min(annCorr(:,1:2), [], 2) - max(annCorr(:,3:end), [], 2));
if isempty(gridness) || isnan(gridness); gridness = NaN; mask = false(sz); rOutBest = NaN; return; end
rOutBest = radii(iBest);
mask = ~cMask & distMap < rOutBest;
end
