function G = gridPeaksFromTroughs(sac, options)
% gridPeaksFromTroughs - find the 6 inner peaks of a grid SAC via the troughs
% around the central peak.
%
%       G = scanpix.analysis.gridPeaksFromTroughs(sac);
%
% package: scanpix.analysis (used by scanpix.analysis.gridprops_v2)
%
% Steps
%  1) polar resampling of the SAC around its centre.
%  2) extent of central peak r0 = first minimum of angle averaged profile;
%     initial axis direction phi from the 6-fold Fourier phase in r0..2.5*r0.
%  3) for each of the 3 independent axes (SAC is point symmetric, axis k+3 = axis k):
%       a) trough along the axis = deepest point before the profile rises again by
%          'prom' (so small shoulders on the flank of the central peak are skipped)
%       b) peak = maximum of the SAC in a sector (+-'sector' deg) beyond the trough.
%          If the maximum sits at the NaN edge of the SAC (profile still rising)
%          the peak is flagged 'censored' (true peak at or beyond the edge).
%       c) axis direction is set to the direction of the (uncensored) peak, and
%          a)-b) are repeated ('nIter' times), so axes needn't be 60 deg apart.
%  4) scale = 2 * mean trough distance (trough sits at half the lattice vector
%     along the line centre->peak for ANY lattice, i.e. also elliptical grids)
%  5) ellipses with centre fixed at the origin (A x^2 + B xy + C y^2 = 1):
%       - ellPeak:   through the 6 peaks (censored peaks replaced by 2 * trough)
%       - ellTrough: through 2 * the 6 trough vectors
%
% Output G (all distances in SAC bins, x,y in image coords, y down)
%   .scale            2 * mean trough distance
%   .scalePeak        mean distance of uncensored peaks
%   .orientation      gridprops-style orientation (rad) from peak axes
%   .peakXY [6x2]     peak coordinates; .peakR, .peakAng, .peakVal, .censored [6x1]
%   .troughXY [6x2]   trough coordinates; .troughR, .troughVal, .troughOK (rise found)
%   .axisAngle [6x1]  final axis directions
%   .ellPeak, .ellTrough  struct a (major semi axis), b (minor), orient (rad), ratio a/b
%   .nCensored
%
% LM / Claude 2026

arguments
    sac {mustBeNumeric}
    options.prom (1,1) {mustBeNumeric} = 0.1;     % rise (r) needed to confirm a trough
    options.sector (1,1) {mustBeNumeric} = 20;    % half width (deg) of peak search sector
    options.maxDev (1,1) {mustBeNumeric} = 20;    % max deviation (deg) of an axis from regular hex
    options.wedge (1,1) {mustBeNumeric} = 3;      % angular average (deg) for trough profiles
    options.nIter (1,1) {mustBeNumeric} = 2;
    options.minR (1,1) {mustBeNumeric} = 2;
    options.edgeTol (1,1) {mustBeNumeric} = 1;    % bins from NaN edge counting as 'at edge'
    options.peakDrop (1,1) {mustBeNumeric} = 0.05; % a peak must be followed by this drop (r) before the edge
    options.troughWin (1,2) {mustBeNumeric} = [0.6 1.4]; % axis troughs searched in this range x r0
    options.band (1,2) {mustBeNumeric} = [1.2 3.2];      % candidate peaks searched in this range x r0
    options.cRange (1,2) {mustBeNumeric} = [0.75 1.2];   % allowed |P| / (2*trough) for a consistent peak (edge peaks read short)
    options.minSep (1,1) {mustBeNumeric} = 35;           % min angle (deg) between two lattice axes
    options.minProm (1,1) {mustBeNumeric} = 0.1;         % a candidate peak must exceed its ray trough by this (r)
    options.minContrast (1,1) {mustBeNumeric} = -inf;    % reject lattice if peak-vs-hole contrast is below this
    options.fallback (1,1) {mustBeNumericOrLogical} = true;  % hex-template fallback if no lattice is accepted
    options.require3 (1,1) {mustBeNumericOrLogical} = false; % reject lattice if the third axis has no detected peak
end

nanEll = struct('a', NaN, 'b', NaN, 'orient', NaN, 'ratio', NaN);
G = struct('scale', NaN, 'scalePeak', NaN, 'orientation', NaN, 'peakXY', nan(6,2), 'peakR', nan(6,1), ...
    'peakAng', nan(6,1), 'peakVal', nan(6,1), 'censored', true(6,1), 'troughXY', nan(6,2), 'troughR', nan(6,1), ...
    'troughVal', nan(6,1), 'troughOK', false(6,1), 'axisAngle', nan(6,1), 'ellPeak', nanEll, 'ellTrough', nanEll, ...
    'nCensored', NaN, 'r0', NaN, 'method', '', 'contrast', NaN, 'r0a', NaN, 'found3', false);
if isempty(sac) || all(isnan(sac(:))); return; end

% --- polar resampling ----------------------------------------------------
dTh = 1; dr = 0.25;
c    = ceil(size(sac)/2);                                   % row, col
rMax = sqrt(sum((size(sac)/2).^2));                         % reach corners
r    = 0:dr:rMax;
th   = (0:dTh:360-dTh) * pi/180;
[R, TH] = meshgrid(r, th);
Pol  = interp2(sac, c(2) + R.*cos(TH), c(1) + R.*sin(TH), 'linear');   % NaN outside

% --- central peak extent + initial orientation ---------------------------
prof = mean(Pol(:, r <= min(size(sac))/2 - 1), 1, 'omitnan');
r0a  = firstMin(prof, r, options.minR);                    % first dip (can be a shoulder)
r0   = promTrough(prof, r, options.minR, 0.05);            % robust to shoulders, but can skip a weak first ring
if isnan(r0); r0 = r0a; end
if isnan(r0a); r0a = r0; end
G.r0a = r0a;
G.r0 = r0;
if isnan(r0); return; end

% --- primary: lattice search on trough-consistent peaks ------------------
L = latticeSearch(sac, Pol, th, r, [min(r0a, r0) max(r0a, r0)], c, options);
G.method = 'lattice';
if isfield(L, 'cand'); G.cand = L.cand; end
G.contrast = L.contrast; G.found3 = L.found3;
if ~L.ok && ~options.fallback
    G.method = 'none'; return
end
if L.ok
    ax = L.pA; tR = L.tR; tV = L.tV; tOK = L.tOK; pR = L.pR; pA = L.pA; pV = L.pV; cens = L.cens;
else
G.method = 'hexTemplate';
annR = r >= r0 & r <= min(2.5*r0, min(size(sac))/2 - 1);
A    = mean(Pol(:, annR), 2, 'omitnan');  A(isnan(A)) = mean(A, 'omitnan');
phi  = angle(sum((A - mean(A)) .* exp(1i*6*th(:)))) / 6;

ax0  = phi + (0:2)'*pi/3;              % regular-hex prediction for the 3 axes
% point each axis at the strongest SAC values in the band where the first-ring peaks are
% expected (peak at ~2*axis trough ~ 1.9*r0), within maxDev of the hex prediction
band = r >= 1.3*r0 & r <= 2.6*r0;
B    = max(Pol(:, band), [], 2, 'omitnan');
B    = mean([B, circshift(B, -round(180/dTh))], 2, 'omitnan');          % point symmetry
B(isnan(B)) = -inf;
Bs   = movmean([B(end-2:end); B; B(1:3)], 5); Bs = Bs(4:end-3);          % 5 deg circular smoothing
ax   = ax0;
for k = 1:3
    dA = angle(exp(1i*(th(:) - ax0(k))));
    w  = find(abs(dA) <= options.maxDev*pi/180);
    if any(isfinite(Bs(w)))
        [~, j] = max(Bs(w)); ax(k) = ax0(k) + dA(w(j));
    end
end
ax0  = ax;                             % deviations in the refinement loop are relative to this
[tR, tV, tOK, pR, pA, pV, cens] = deal(nan(3,1));
for it = 1:options.nIter+1
    for k = 1:3
        % a) trough along axis
        p = wedgeProfile(Pol, th, ax(k), options.wedge);
        [tR(k), tV(k), tOK(k)] = axisTrough(p, r, r0, options);
        if isnan(tR(k)); continue; end
        % b) peak in sector beyond trough (search both half planes; SAC point symmetric -> average)
        [pR(k), pA(k), pV(k), cens(k)] = sectorPeak(Pol, th, r, ax(k), tR(k), tV(k), options);
    end
    if it > options.nIter; break; end
    % c) update axis directions from uncensored peaks (limit deviation from hex prediction)
    for k = 1:3
        if ~isnan(pA(k)) && ~cens(k)
            d = angle(exp(1i*(pA(k) - ax0(k))));
            if abs(d) <= options.maxDev*pi/180; ax(k) = ax0(k) + d; end
        end
    end
end
end   % hex-template fallback

% --- expand to 6 (point symmetry) and pack -------------------------------
ax6  = mod([ax; ax+pi], 2*pi);
dup  = @(v) [v; v];
G.axisAngle = ax6;
G.troughR   = dup(tR);  G.troughVal = dup(tV);  G.troughOK = dup(tOK == 1);
G.troughXY  = [c(2) + G.troughR.*cos(ax6), c(1) + G.troughR.*sin(ax6)];
pA6 = mod([pA; pA+pi], 2*pi);
G.peakR = dup(pR); G.peakAng = pA6; G.peakVal = dup(pV); G.censored = dup(cens == 1);
G.peakXY  = [c(2) + G.peakR.*cos(pA6), c(1) + G.peakR.*sin(pA6)];
G.nCensored = sum(cens == 1);
if sum(~isnan(tR)) >= 2; G.scale = 2*mean(tR, 'omitnan'); end
if any(cens == 0);       G.scalePeak = mean(pR(cens == 0)); end

% orientation, gridprops style (3 peaks with smallest |angle|, y flipped)
if all(~isnan(pA6))
    o = atan2(sin(pA6), cos(pA6));
    [~, s] = sort(abs(o));
    G.orientation = circ_mean_(-o(s(1:3)));
end

% --- ellipses (centre fixed) ---------------------------------------------
v  = [pR.*cos(pA), pR.*sin(pA)];
vt = 2*[tR.*cos(ax), tR.*sin(ax)];
v(cens == 1 | isnan(pR), :) = vt(cens == 1 | isnan(pR), :);   % censored -> 2*trough
G.ellPeak   = fitCentredEllipse(v);
G.ellTrough = fitCentredEllipse(vt);
end


% =========================================================================
function L = latticeSearch(sac, Pol, th, r, r0, c, o)
L.contrast = NaN; L.found3 = false;
% Candidate peaks = regional maxima of the SAC in band*r0. A candidate P is consistent
% if the deepest minimum along the ray centre->P sits at |P|/2 (|log(|P|/2tr)| <= maxLogC)
% and no other peak lies on the ray before it. Then pick the pair (a1, a2) + third vector
% a3 = shorter of a1+-a2 that best forms a lattice (consistency, a3 present, small radius).
L.ok = false;
S = sac; nanM = isnan(S); S(nanM) = -inf;
BW = imregionalmax(S) & ~nanM;
CC = bwconncomp(BW, 8);
if CC.NumObjects < 2; return; end
nearNaN = imdilate(nanM, strel('square', 2*o.edgeTol + 1));
[rows, cols] = deal(nan(CC.NumObjects,1)); [val, cen] = deal(nan(CC.NumObjects,1));
for k = 1:CC.NumObjects
    [ri, ci] = ind2sub(size(S), CC.PixelIdxList{k});
    [val(k), j] = max(sac(CC.PixelIdxList{k}));
    rr = ri(j); cc = ci(j);
    cen(k) = any(nearNaN(CC.PixelIdxList{k}));
    if ~cen(k) && rr > 1 && cc > 1 && rr < size(S,1) && cc < size(S,2)
        w = sac(rr-1:rr+1, cc-1:cc+1); w = w - min(w(:)); [X, Y] = meshgrid(cc-1:cc+1, rr-1:rr+1);
        rows(k) = sum(w(:).*Y(:))/sum(w(:)); cols(k) = sum(w(:).*X(:))/sum(w(:));
    else
        rows(k) = rr; cols(k) = cc;
    end
end
v   = [cols - c(2), rows - c(1)];
rad = sqrt(sum(v.^2, 2));
in  = rad >= o.band(1)*r0(1) & rad <= o.band(2)*r0(end);   % r0 = [first dip, shoulder-proof trough]
v = v(in,:); rad = rad(in); val = val(in); cen = cen(in);
if numel(rad) < 2; return; end
% consistency of each candidate with the trough along its ray
n = numel(rad); [tr, tv, logC, okC] = deal(nan(n,1));
for k = 1:n
    a = atan2(v(k,2), v(k,1));
    p = wedgeProfile(Pol, th, a, o.wedge);
    w = find(r >= 0.3*rad(k) & r <= 0.8*rad(k));
    w = w(w > 1 & w < numel(p));
    w = w(~isnan(p(w-1)) & ~isnan(p(w+1)) & p(w) <= p(w-1) & p(w) < p(w+1));
    if isempty(w); okC(k) = 0; continue; end
    [tv(k), j] = min(p(w)); tr(k) = refine(p, r, w(j), -1);
    logC(k) = log(rad(k) / (2*tr(k)));
    % no other peak on the ray between trough and P
    seg = p(r > tr(k) & r < rad(k) - 1.5); clean = true;
    for i = 2:numel(seg)-1
        if seg(i) >= seg(i-1) && seg(i) > seg(i+1) && min(seg(i+1:end)) <= seg(i) - o.peakDrop; clean = false; break; end
    end
    okC(k) = logC(k) >= log(o.cRange(1)) && logC(k) <= log(o.cRange(2)) && clean && val(k) > tv(k) + o.minProm;
end
acc = find(okC == 1);
F = griddedInterpolant({1:size(sac,1), 1:size(sac,2)}, sac, 'linear', 'none');   % built once for the lattice scoring
L.cand = table(v(:,1), v(:,2), rad, val, cen, tr, logC, okC, 'VariableNames', {'x','y','r','val','cens','tr','logC','ok'});
% half plane representatives of accepted candidates (point symmetry)
ang = mod(atan2(v(:,2), v(:,1)), pi);
best = -inf; sel = [];
% pair search over one copy of each mirror pair (P and -P give the same lattice hypotheses)
accP = acc;
if numel(acc) > 1
    dm = sqrt((v(acc,1) + v(acc,1)').^2 + (v(acc,2) + v(acc,2)').^2);   % |v_i + v_j|: ~0 for mirror copies
    dm(logical(eye(numel(acc)))) = inf;
    isMirrorLater = any(triu(dm < 1.5, 1), 1);                           % drop the later copy of each pair
    accP = acc(~isMirrorLater);
end
for ii = 1:numel(accP)
    for jj = ii+1:numel(accP)
        i = accP(ii); j = accP(jj);
        sep = abs(angle(exp(2i*(ang(i) - ang(j))))) / 2 * 180/pi;      % 0..90
        if sep < o.minSep; continue; end
        vi = v(i,:); vj = v(j,:);
        if dot(vi, vj) < 0; vj = -vj; end                              % angle < 90
        % third vector vi - vj or vi + vj; snap to a candidate if one is close
        for p3c = {vi - vj, vi + vj}
            p3 = p3c{1};
            if norm(p3) > 1.6*max(norm(vi), norm(vj)); continue; end
            d  = [sqrt(sum((v - p3).^2, 2)), sqrt(sum((v + p3).^2, 2))];
            [d3, li] = min(d(:)); [k, sg] = ind2sub(size(d), li);
            found = d3 <= 0.25*norm(p3);
            if found; p3 = v(k,:) * (3 - 2*sg); end                     % sg 1: +v, 2: -v
            % nearest-neighbour parsimony: consistent peaks closer than this lattice's own
            % peaks that the lattice does not explain (point-symmetric pairs counted once)
            Vh = [vi; vj; p3]; rLat = max(sqrt(sum(Vh.^2, 2)));
            Nn = [Vh; -Vh];
            dn = min(sqrt((v(acc,1) - Nn(:,1)').^2 + (v(acc,2) - Nn(:,2)').^2), [], 2);
            nUn = sum(rad(acc) < 0.9*rLat & dn > 0.25*rad(acc)) / 2;
            sc = latticeContrast(F, c, Vh) - 0.25*nUn;
            if sc > best; best = sc; sel = [i j k found]; p3best = p3; vjbest = vj; end
        end
    end
end
L.contrast = best; L.found3 = ~isempty(sel) && sel(4);
if isempty(sel) || best < o.minContrast || (o.require3 && ~sel(4)); return; end
i = sel(1); j = sel(2); k = sel(3);
V = [v(i,:); vjbest; p3best]; isCen = [cen(i); cen(j); 1]; pv = [val(i); val(j); NaN];
if sel(4); isCen(3) = cen(k); pv(3) = val(k); end
L.pA = atan2(V(:,2), V(:,1)); L.pR = sqrt(sum(V.^2, 2)); L.pV = pv; L.cens = double(isCen);
% troughs along the final rays
[L.tR, L.tV, L.tOK] = deal(nan(3,1));
for q = 1:3
    p = wedgeProfile(Pol, th, L.pA(q), o.wedge);
    w = find(r >= 0.3*L.pR(q) & r <= 0.8*L.pR(q)); w = w(w > 1 & w < numel(p));
    w = w(~isnan(p(w-1)) & ~isnan(p(w+1)) & p(w) <= p(w-1) & p(w) < p(w+1));
    if ~isempty(w); [L.tV(q), jj] = min(p(w)); L.tR(q) = refine(p, r, w(jj), -1); L.tOK(q) = 1;
    else; [L.tR(q), L.tV(q), L.tOK(q)] = axisTrough(p, r, L.pR(q)/2, o); end
end
[~, s] = sort(mod(L.pA, pi));
for f = {'pA','pR','pV','cens','tR','tV','tOK'}; L.(f{1}) = L.(f{1})(s); end
L.ok = true;
end

function sc = latticeContrast(F, c, V)
% mean SAC at the 6 predicted neighbour peaks minus mean SAC at the 6 holes between them
% (centroids of the triangles centre - n_k - n_k+1). NaN (outside SAC) ignored.
N   = [V; -V];
[~, s] = sort(atan2(N(:,2), N(:,1)));  N = N(s,:);
H   = (N + circshift(N, -1)) / 3;
f   = @(P) F(c(1) + P(:,2), c(2) + P(:,1));
pk  = f(N); ho = f(H);
okp = ~isnan(pk); okh = ~isnan(ho);
if sum(okp) < 3 || sum(okh) < 3; sc = -inf; return; end
sc  = sum(pk(okp))/sum(okp) - sum(ho(okh))/sum(okh);
end

% =========================================================================
function p = wedgeProfile(Pol, th, ang, wedge)
dA = abs(angle(exp(1i*(th(:) - ang)))) * 180/pi;
p  = mean(Pol(dA <= wedge, :), 1, 'omitnan');
end

function [rm, vm, ok] = axisTrough(p, r, r0, o)
% deepest local minimum within troughWin * r0 (r0 = trough of the angle averaged profile).
% Picks the real trough over small shoulders on the central peak flank (shallower) and
% never jumps past a weak first ring (the next minimum along the axis is at ~3*r0).
% Falls back to promTrough if there is no local minimum in the window.
rm = NaN; vm = NaN; ok = 0;
inW = find(r >= max(o.minR, o.troughWin(1)*r0) & r <= o.troughWin(2)*r0);
cand = inW(inW > 1 & inW < numel(p));
cand = cand(~isnan(p(cand-1)) & ~isnan(p(cand+1)) & p(cand) <= p(cand-1) & p(cand) < p(cand+1));
if ~isempty(cand)
    [vm, j] = min(p(cand)); rm = refine(p, r, cand(j), -1); ok = 1;
else
    [rm, vm, ok] = promTrough(p, r, o.minR, o.prom);
end
end

function [rm, vm, ok] = promTrough(p, r, minR, prom)
% deepest point before the profile rises 'prom' above the running minimum
rm = NaN; vm = NaN; ok = 0;
i0 = find(r >= minR, 1);
runMin = inf; iMin = NaN;
for i = i0:numel(p)
    if isnan(p(i)); break; end
    if p(i) < runMin; runMin = p(i); iMin = i; end
    if p(i) > runMin + prom; ok = 1; break; end
end
if isnan(iMin) || iMin == i0; return; end
if ~ok && (iMin == numel(p) || isnan(p(iMin+1))); return; end   % still falling at edge
rm = refine(p, r, iMin, -1); vm = p(iMin);
end

function [pr, pa, pv, cens] = sectorPeak(Pol, th, r, ang, tr, tv, o)
pr = NaN; pa = NaN; pv = NaN; cens = 1;
dA  = abs(angle(exp(1i*(th(:) - ang)))) * 180/pi;
dA2 = abs(angle(exp(1i*(th(:) - ang - pi)))) * 180/pi;
iS  = find(dA <= o.sector); iS2 = find(dA2 <= o.sector);
[~, ord] = sort(th(iS)); iS = iS(ord);
[~, ord] = sort(mod(th(iS2) - pi, 2*pi)); iS2 = iS2(ord);
% the SAC is point symmetric: average both half planes (also fills NaNs on one side)
n   = min(numel(iS), numel(iS2));
S   = mean(cat(3, Pol(iS(1:n), :), Pol(iS2(1:n), :)), 3, 'omitnan');
S(:, r <= tr | r > 3*tr) = NaN;     % first ring only (peak expected near 2*tr; 2nd ring at ~4*tr)
if all(isnan(S(:))); return; end
% along each ray: first local maximum beyond the trough. A ray that keeps rising
% until the NaN edge (or the 3*tr limit) only yields an edge candidate (censored).
nA = size(S,1);
[cR, cV, cC] = deal(nan(nA,1));
for a = 1:nA
    v  = find(~isnan(S(a,:)));
    if numel(v) < 3; continue; end
    p  = S(a,:);
    iL = v(end);
    for i = v(2):iL-1
        if p(i) >= p(i-1) && p(i) > p(i+1) && (r(iL) - r(i)) >= o.edgeTol && min(p(i+1:iL)) <= p(i) - o.peakDrop
            cR(a) = i; cV(a) = p(i); cC(a) = 0; break
        end
    end
    if isnan(cR(a)); [cV(a), j] = max(p(v)); cR(a) = v(j); cC(a) = 1; end
end
% prefer proper peaks over edge candidates
if any(cC == 0); cV(cC ~= 0) = -inf; end
if ~any(isfinite(cV)); return; end
% the trough sits at half the lattice vector -> peak expected at ~2*tr. Among candidates that
% are reasonably strong (>= half way from trough to best candidate) take the one closest to 2*tr,
% then snap to the highest candidate within 2 bins of it (rays through the peak centre).
ok   = isfinite(cV);
rC   = nan(nA,1); rC(ok) = r(cR(ok));
aC   = th(iS(1:nA)); aC = aC(:);
best = max(cV(ok));
elig = ok & cV >= tv + 0.5*(best - tv);
[~, j0] = min(abs(rC - 2*tr) + 1e6*~elig);
xy   = [rC.*cos(aC), rC.*sin(aC)];
near = ok & sqrt(sum((xy - xy(j0,:)).^2, 2)) <= 2;
cV(~near) = -inf;
[pv, ia] = max(cV);
ir = cR(ia); cens = cC(ia);
pa = th(iS(ia)); pr = r(ir);
if ~cens; pr = refine(S(ia,:), r, ir, 1); end
end

function rr = refine(p, r, i, sgn)
rr = r(i);
if i > 1 && i < numel(p) && all(~isnan(p(i-1:i+1)))
    den = p(i-1) - 2*p(i) + p(i+1);
    if sgn*den < 0          % min: den > 0, max: den < 0
        rr = r(i) + 0.5*(p(i-1) - p(i+1))/den * (r(2) - r(1));
    end
end
end

function rm = firstMin(p, r, minR)
rm = NaN;
for i = find(r >= minR, 1):numel(p)-1
    if isnan(p(i+1)) || isnan(p(i)); return; end
    if p(i) <= p(i-1) && p(i) < p(i+1); rm = r(i); return; end
end
end

function E = fitCentredEllipse(v)
% A x^2 + B xy + C y^2 = 1 through points v (rows x,y); point symmetry -> centre 0
E = struct('a', NaN, 'b', NaN, 'orient', NaN, 'ratio', NaN);
v = v(all(isfinite(v), 2), :);
if size(v,1) < 3; return; end
M = [v(:,1).^2, v(:,1).*v(:,2), v(:,2).^2];
q = M \ ones(size(v,1), 1);
Q = [q(1) q(2)/2; q(2)/2 q(3)];
[V, D] = eig(Q); d = diag(D);
if any(d <= 0); return; end
[d, s] = sort(d); V = V(:, s);              % smallest eigenvalue -> major axis
E.a = 1/sqrt(d(1)); E.b = 1/sqrt(d(2));
E.orient = atan2(V(2,1), V(1,1));
E.ratio  = E.a / E.b;
end

function m = circ_mean_(a)
m = angle(mean(exp(1i*a)));
end
