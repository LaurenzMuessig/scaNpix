function exportViaPython( hCanvas, fNameOut, resolution, pythonExe )
% Export a figure/uipanel of mixed axes (images, line plots, area plots,
% text - e.g. a scanpix.plot.multPlot grid mixing plotRateMap/plotDirMap/
% plotSpeedMap/plotLinMaps panels) to a genuinely vector PDF by handing the
% panel data to matplotlib.
% package: scanpix.helpers
%
% MATLAB's own vector PDF export (exportgraphics/print) currently has no
% reliable way to embed rasterised content (imagesc heatmaps etc) correctly -
% see the 'contentType' options in scanpix.helpers.saveFigAsPDF. matplotlib's
% PDF backend does not share that bug and gives real vector text, so this
% function: (1) walks every axes and every child graphics object in it,
% converts what it recognises (image/line/area/text) into a plain struct
% plus the axis-level formatting (limits, ticks, labels, direction,
% aspect), (2) saves that to a .mat file, (3) shells out to
% plotMultiFromMat.py, which rebuilds the same layout in matplotlib and
% writes the PDF.
%
% Not supported (silently skipped, with a warning): polaraxes/polarplot,
% surfaces, patches, legends, colorbars, scatter, bar.
%
% Usage:    scanpix.helpers.exportViaPython( hCanvas, fNameOut, resolution, pythonExe )
%
% Inputs:   hCanvas    - handle to the uipanel or figure containing the axes
%           fNameOut   - full output .pdf path
%           resolution - DPI for any embedded raster (image) panels
%           pythonExe  - python executable (must have scipy + matplotlib)
%
% LM 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hCanvas.Units = 'pixels';
canvasPos = hCanvas.Position; % [x y w h], defines the output page size

% hCanvas.Children (like any Children list) is stored newest-first; flip to
% get creation order, so stacked/overlaid sibling axes (e.g.
% scanpix.plot.plotGridProps's transparent foreground axes on top of a
% background axes at the same position) are drawn back-to-front correctly
allKids  = flipud(hCanvas.Children);
hAx      = allKids(arrayfun(@(h) isa(h,'matlab.graphics.axis.Axes'), allKids));
axesData = struct('left',{},'bottom',{},'width',{},'height',{}, ...
    'xlim',{},'ylim',{},'xdir',{},'ydir',{},'visible',{},'aspectSquare',{},'facecolor',{}, ...
    'xlabel',{},'ylabel',{},'title',{},'xtick',{},'ytick',{},'xticklabel',{},'yticklabel',{}, ...
    'children',{});
nAx = 0;

for k = 1:numel(hAx)
    ax = hAx(k);
    ax.Units = 'pixels';
    axPos    = ax.Position; % [x y w h], relative to hCanvas, pixels

    children = {};
    hKids = flipud(ax.Children); % Children is back-to-front reversed; flip to get draw order
    for j = 1:numel(hKids)
        h = hKids(j);
        switch h.Type
            case 'image'
                cdata   = double(h.CData);
                nanMask = isnan(cdata); % MATLAB renders NaN CData as transparent regardless of AlphaData - must replicate that, not just clamp the colour index
                cmap    = colormap(ax);
                climAx  = clim(ax);
                nColors = size(cmap,1);
                % replicate MATLAB's own scaled-CData-to-colormap mapping
                % (see 'help image' CDataMapping) so the RGB matplotlib
                % gets matches what MATLAB itself shows
                idxImg = floor( (cdata - climAx(1)) / (climAx(2)-climAx(1)) * nColors ) + 1;
                idxImg = max(1, min(nColors, idxImg));
                idxImg(nanMask) = 1; % arbitrary valid index - alpha below is forced to 0 here anyway

                alphaRaw = double(h.AlphaData);
                if isscalar(alphaRaw)
                    alphaFull = alphaRaw * ones(size(cdata));
                else
                    alphaFull = alphaRaw;
                end
                alphaFull(nanMask) = 0;

                c = struct('type','image','rgb',reshape(cmap(idxImg(:),:), [size(idxImg) 3]), ...
                    'alpha',alphaFull, ...
                    'xdata',double(h.XData([1 end])),'ydata',double(h.YData([1 end])));
            case 'line'
                c = struct('type','line','xdata',double(h.XData),'ydata',double(h.YData), ...
                    'color',h.Color,'linewidth',h.LineWidth,'linestyle',h.LineStyle);
            case 'area'
                c = struct('type','area','xdata',double(h.XData),'ydata',double(h.YData), ...
                    'facecolor',h.FaceColor,'edgecolor',h.EdgeColor,'basevalue',h.BaseValue);
            case 'text'
                textUnits = h.Units;
                pos = h.Position;
                if strcmp(textUnits,'data')
                    c = struct('type','text','mode','data','x',pos(1),'y',pos(2),'string',{h.String}, ...
                        'color',h.Color,'fontsize',h.FontSize,'ha',h.HorizontalAlignment,'va',h.VerticalAlignment);
                else
                    % 'pixels'/'normalized'/etc - convert to an absolute figure-fraction
                    % position (matching how axPos itself is converted below) rather than
                    % trying to replicate every MATLAB text unit in matplotlib
                    if strcmp(textUnits,'pixels')
                        absX = axPos(1) + pos(1); absY = axPos(2) + pos(2);
                    else % 'normalized' or similar - treat as a fraction of the axes box
                        absX = axPos(1) + pos(1)*axPos(3); absY = axPos(2) + pos(2)*axPos(4);
                    end
                    c = struct('type','figtext','x',absX/canvasPos(3),'y',absY/canvasPos(4),'string',{h.String}, ...
                        'color',h.Color,'fontsize',h.FontSize,'ha',h.HorizontalAlignment,'va',h.VerticalAlignment);
                end
            otherwise
                warning('scanpix:helpers:exportViaPython:unsupportedType', ...
                    '''%s'' objects are not supported by the python export path and will be skipped.', h.Type);
                continue
        end
        children{end+1} = c; %#ok<AGROW>
    end

    if isempty(children)
        continue % nothing recognised in this axes - skip it entirely
    end

    nAx = nAx+1;
    axesData(nAx).left   = axPos(1) / canvasPos(3);
    axesData(nAx).bottom = axPos(2) / canvasPos(4);
    axesData(nAx).width  = axPos(3) / canvasPos(3);
    axesData(nAx).height = axPos(4) / canvasPos(4);
    axesData(nAx).xlim   = ax.XLim;
    axesData(nAx).ylim   = ax.YLim;
    axesData(nAx).xdir   = ax.XDir;
    axesData(nAx).ydir   = ax.YDir;
    axesData(nAx).visible = char(ax.Visible); % ax.Visible is a matlab.lang.OnOffSwitchState object (not plain char) in current MATLAB - save() can't serialise it usefully, so it arrives in Python as an unusable opaque blob unless forced to char here
    axesData(nAx).facecolor    = ax.Color; % 'none' for a transparent overlay axes (e.g. plotGridProps) - must not paint over whatever is stacked beneath it
    axesData(nAx).aspectSquare = isequal(ax.DataAspectRatio,[1 1 1]) || isequal(ax.PlotBoxAspectRatio,[1 1 1]);
    axesData(nAx).xlabel  = get(get(ax,'xlabel'),'String');
    axesData(nAx).ylabel  = get(get(ax,'ylabel'),'String');
    axesData(nAx).title   = get(get(ax,'title'),'String');
    axesData(nAx).xtick   = ax.XTick;
    axesData(nAx).ytick   = ax.YTick;
    axesData(nAx).xticklabel = ax.XTickLabel;
    axesData(nAx).yticklabel = ax.YTickLabel;
    axesData(nAx).children   = children;
end

if nAx == 0
    error('scanpix:helpers:exportViaPython:noSupportedContent','No axes with supported content (image/line/area/text) were found - nothing to export.');
end

figW = canvasPos(3)/96; % inches, treating the panel's pixel layout as a 96 dpi source
figH = canvasPos(4)/96;

matFile = [tempname '.mat'];
save(matFile,'axesData','figW','figH','-v7');
cleanupMat = onCleanup(@() delete(matFile));

scriptPath = fullfile(fileparts(mfilename('fullpath')),'plotMultiFromMat.py');
cmd = sprintf('"%s" "%s" "%s" "%s" %g', pythonExe, scriptPath, matFile, fNameOut, resolution);
[status, cmdOut] = system(cmd);
if status ~= 0
    error('scanpix:helpers:exportViaPython:pythonFailed', ...
        'Python export failed (exit %d). Check pythonExe points to a Python with scipy+matplotlib installed.\n%s', status, cmdOut);
end

end
