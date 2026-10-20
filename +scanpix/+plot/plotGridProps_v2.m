function Props = plotGridProps_v2(autoCorr, options)
% plotGridProps_v2 - plot grid properties from scanpix.analysis.gridprops_v2
% package: scanpix.plot
%
%       scanpix.plot.plotGridProps_v2(autoCorr);
%       Props = scanpix.plot.plotGridProps_v2(autoCorr, 'binSize', 2.5);
%       prms  = scanpix.maps.defaultParamsRateMaps;
%       scanpix.plot.plotGridProps_v2(autoCorr, 'binSize', prms.rate.binSizeSpat, 'gpParams', prms.gridPropsV2);
%
% Two panels: (1) SAC, (2) ellipse corrected SAC. Only the annulus used for the
% gridness score is shown in colour (rest in grey). Overlays:
%   cyan dots        - troughs between central peak and the six peaks
%   white circles    - peaks used (squares = peak cut off by the SAC edge, distance from 2x trough)
%   black lines      - the 3 grid axes (peaks used for orientation)
%   white dashed     - cartesian reference; red arc = grid offset
%   magenta dashed   - ellipse through the six peaks (panel 1)
% All properties are listed next to each panel.
%
% Inputs
%   autoCorr          - spatial autocorrelogram
%   'axArray'         - (1,2) cell of axes (default: new scanpix.plot.multPlot figure)
%   'binSize'         - bin size in cm; if given, distances are also shown in cm
%   'gpParams'        - struct of name-value options for scanpix.analysis.gridprops_v2
%                       (e.g. prms.gridPropsV2 from scanpix.maps.defaultParamsRateMaps)
%   'Props'           - precomputed output of gridprops_v2 (skips the computation)
%
% Output
%   Props             - gridprops_v2 output structure
%
% LM / Claude 2026

arguments
    autoCorr {mustBeNumeric}
    options.axArray cell = {}
    options.binSize {mustBeScalarOrEmpty} = []
    options.gpParams (1,1) struct = struct()
    options.Props (1,1) struct = struct()
end

%% compute
if isempty(fieldnames(options.Props))
    gpArgs = namedargs2cell(options.gpParams);
    [~, Props] = scanpix.analysis.gridprops_v2(autoCorr, gpArgs{:});
else
    Props = options.Props;
end
if isempty(options.axArray)
    options.axArray = scanpix.plot.multPlot([1 2], 'plotsize', [260 260], 'plotsep', [230 40], 'figname', 'gridprops_v2');
end

%% panels
% properties per panel (panel 2 = ellipse corrected SAC)
Q(1) = struct('gridness', Props.gridness, 'wavelength', Props.wavelength, 'orientation', Props.orientation, ...
    'orientationFull', Props.orientationFull, 'offset', Props.offset, 'offsetFull', Props.offsetFull, ...
    'closestPeaksCoord', Props.closestPeaksCoord, 'pk', Props.peakCoords(2:end,:), 'tr', Props.troughCoords);
Q(2) = struct('gridness', Props.gridnessEll, 'wavelength', Props.wavelengthEll, 'orientation', Props.orientationEll, ...
    'orientationFull', Props.orientationFullEll, 'offset', Props.offsetEll, 'offsetFull', Props.offsetFullEll, ...
    'closestPeaksCoord', Props.closestPeaksCoordEll, 'pk', Props.peakCoordsEll, 'tr', Props.troughCoordsEll);
sacs   = {autoCorr, Props.acReg};
masks  = {Props.gridMask, Props.gridMaskEll};
titles = {'spatial autocorr', 'spatial autocorr - ellipse corrected'};
for p = 1:2
    ax  = options.axArray{p};
    sac = sacs{p};
    cla(ax);
    if isempty(sac) || all(isnan(sac(:)))
        axis(ax, 'off'); title(ax, [titles{p} ' (n/a)']); continue
    end
    sz = size(sac); c = ceil(sz/2);

    % a. annulus in colour (jet), rest in grey - composed as one RGB image
    lims = [min(sac(:), [], 'omitnan') 1];
    if ~(lims(1) < lims(2)); lims = [-1 1]; end
    idx  = @(x, n) max(1, min(n, round(1 + (x - lims(1)) ./ diff(lims) * (n - 1))));
    nC   = 256; jc = jet(nC); gc = gray(nC);
    I    = idx(sac, nC); I(isnan(sac)) = 1;
    rgb  = reshape(gc(I(:), :), [sz 3]);
    if any(masks{p}(:))
        cj = reshape(jc(I(:), :), [sz 3]);
        m3 = repmat(masks{p}, [1 1 3]);
        rgb(m3) = cj(m3);
    end
    nanM = repmat(isnan(sac), [1 1 3]); rgb(nanM) = 0;      % NaN bins (outside SAC / low overlap) = black
    image(ax, rgb); axis(ax, 'image', 'off'); hold(ax, 'on');

    % b. annulus edges
    th = linspace(0, 2*pi, 180);
    ann = Props.annulus; if p == 2; ann = Props.annulusEll; end
    for r = ann(isfinite(ann))
        plot(ax, c(2) + r*cos(th), c(1) + r*sin(th), '-', 'Color', [1 1 1 0.8], 'LineWidth', 1);
    end

    % c. cartesian reference
    plot(ax, [0.5 sz(2)+0.5], c([1 1]), '--w', 'LineWidth', 1.5);
    plot(ax, c([2 2]), [0.5 sz(1)+0.5], '--w', 'LineWidth', 1.5);

    if Props.latticeFound
        q = Q(p); pk = q.pk; tr = q.tr;
        cens = [Props.censored; Props.censored] == 1;
        % d. grid axes (3 peaks with smallest |angle|, as used for orientation)
        if all(isfinite(q.closestPeaksCoord(:)))
            for k = 1:3; plot(ax, [c(2) q.closestPeaksCoord(k,1)], [c(1) q.closestPeaksCoord(k,2)], 'k', 'LineWidth', 2); end
        end
        % e. troughs + peaks
        if any(isfinite(tr(:)))
            plot(ax, tr(:,1), tr(:,2), 'o', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
        end
        if any(isfinite(pk(:)))
            plot(ax, pk(~cens,1), pk(~cens,2), 'o', 'Color', 'w', 'MarkerSize', 11, 'LineWidth', 2.5);
            plot(ax, pk(cens,1),  pk(cens,2),  's', 'Color', 'w', 'MarkerSize', 11, 'LineWidth', 2);
        end
        % f. ellipse through the peaks (panel 1) / the circle it is mapped to (panel 2)
        if all(isfinite(Props.ellAbScale))
            a = Props.ellAbScale(1); b = Props.ellAbScale(2); o = Props.ellOrient;
            if p == 2; a = b; end
            xy = [cos(o) -sin(o); sin(o) cos(o)] * [a*cos(th); b*sin(th)];
            plot(ax, c(2) + xy(1,:), c(1) + xy(2,:), 'm--', 'LineWidth', 1.2);
        end
        % g. offset arc (same convention as scanpix.plot.plotGridProps)
        if all(isfinite(q.offsetFull))
            [minOffset, ind] = min(q.offsetFull);
            if ind == 1
                if q.orientationFull(ind) > 0; shift = 0; else; shift = 2*pi - minOffset; end
            else
                if q.orientationFull(ind) > 0; shift = pi/2 - minOffset; else; shift = 3/2*pi; end
            end
            tha = shift:0.05:minOffset + shift;
            [x, y] = pol2cart(tha, q.wavelength * 0.6);
            plot(ax, x + c(2), c(1) - y, '-r', 'LineWidth', 3);
        end
    end
    xlim(ax, [0.5 sz(2)+0.5]); ylim(ax, [0.5 sz(1)+0.5]);
    hold(ax, 'off');
    title(ax, titles{p});

    % h. property text next to the panel
    text(ax, sz(2) + 2, 1, propText(Props, Q(p), p, options.binSize), 'VerticalAlignment', 'top', ...
        'FontSize', 8.5, 'Interpreter', 'tex', 'Clipping', 'off');
end
end

% -------------------------------------------------------------------------
function s = propText(P, q, p, binSize)
f = @(v, d) num2str(round(v, d));
dist = @(v) local_dist(v, binSize);
L = {['\bfgridness: ' f(q.gridness, 3) '\rm'], ''};
if P.latticeFound && all(isfinite([q.wavelength q.orientation q.offset]))
    L = [L, {sprintf('scale: %s', dist(q.wavelength)), ...
             sprintf('orientation: %s%s', f(rad2deg(q.orientation), 1), char(176)), ...
             sprintf('offset: %s%s', f(rad2deg(q.offset), 1), char(176))}];
else
    L = [L, {'scale: n/a', 'orientation: n/a', 'offset: n/a'}];
end
if p == 2 && P.latticeFound && ~P.ellCorrected
    L = [L, {'', sprintf('a/b %s, weak lattice:', f(P.ellRatio, 2)), 'correction not used,', ...
             'gridness = uncorrected', sprintf('(stretched: %s)', f(P.gridnessEllStretched, 2))}];
end
s = strjoin(L, newline);
end

% -------------------------------------------------------------------------
function s = local_dist(v, binSize)
if isempty(binSize)
    s = sprintf('%.1f bins', v);
else
    s = sprintf('%.1f cm', v * binSize);
end
end
