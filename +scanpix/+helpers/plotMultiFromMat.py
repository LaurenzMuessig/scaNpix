"""
Rebuild a MATLAB axes grid (images, line plots, area plots, text - e.g. a
mixed scanpix.plot.multPlot grid) in matplotlib and save as a PDF.

Exists to work around a MathWorks bug (#3668256, unresolved as of R2026a) where
exportgraphics/print corrupt or blur rasterised content (e.g. imagesc heatmaps)
embedded in vector PDF/EPS output. matplotlib's PDF backend does not share this
bug: raster content is embedded at a resolution set directly by --dpi, and
lines/text stay genuine vector objects.

Called from scanpix.helpers.exportViaPython.m (in turn from
scanpix.helpers.saveFigAsPDF(...,'contentType','python')).

Usage: python plotMultiFromMat.py <matFile> <outPdf> <dpi>
"""
import sys
import scipy.io
import matplotlib
matplotlib.use('Agg')
# matplotlib's PDF backend normally zlib-compresses embedded raster images
# using a PNG-style predictor filter (/DecodeParms << /Predictor 15 >>). That
# renders correctly in most PDF readers (confirmed here: MuPDF/PyMuPDF and
# Claude's own PDF reader both display the file's actual embedded colours
# correctly) but Adobe Illustrator has shown a real, reproducible mismatch on
# these files - wrong, image-specific colours, not a uniform tint - strongly
# suggesting a predictor-decoding incompatibility rather than a data bug.
# Disabling compression removes the predictor entirely (larger files, but
# universally safe raw image data).
matplotlib.rcParams['pdf.compression'] = 0
# matplotlib defaults to embedding text as Type 3 fonts (pdf.fonttype 3) -
# each glyph is its own tiny standalone drawing program rather than part of a
# real font. Illustrator can't treat that as continuous text and explodes it
# into one object per letter. Type 42 embeds the actual TrueType font program,
# which Illustrator (and other PDF editors) read back as normal, single,
# editable text objects.
matplotlib.rcParams['pdf.fonttype'] = 42
# Use Arial instead of matplotlib's default (DejaVu Sans), which Illustrator
# doesn't ship with - Arial is a standard system font so this should always
# be found by matplotlib's font manager.
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
import matplotlib.pyplot as plt
import numpy as np


def as_str(val):
    """MATLAB char/cellstr round-trips through scipy as str, or a 0/N-element ndarray of str.
    A text object whose String was set with an embedded '\\n' comes back from MATLAB as
    MULTIPLE separate per-line strings (matlab.graphics.primitive.Text splits it internally),
    not one string containing a literal newline - rejoin with '\\n' rather than keeping only
    the first line."""
    if isinstance(val, np.ndarray):
        if val.size == 0:
            return ''
        if val.size > 1:
            return '\n'.join(str(v) for v in val.tolist())
        val = val.flat[0]
    return str(val) if val is not None else ''


def as_str_list(val):
    if isinstance(val, np.ndarray):
        return [as_str(v) for v in val.tolist()] if val.size else []
    return [as_str(val)]


def mpl_va(matlab_va):
    """MATLAB VerticalAlignment values don't all match matplotlib's va values."""
    return {'middle': 'center', 'cap': 'top'}.get(matlab_va, matlab_va)


def img_extent(xdata, ydata, shape):
    """MATLAB image XData/YData give the data-coordinate centres of the first/last
    pixel; matplotlib's imshow extent wants the outer edges. Using the real MATLAB
    coordinates (rather than matplotlib's default 0-based pixel-index extent) keeps
    any lines/text drawn in that same data space (e.g. plotGridProps) in registration
    with the image, not off by the 1-based/0-based indexing difference."""
    nrows, ncols = shape[0], shape[1]
    dx = (xdata[1]-xdata[0])/(ncols-1) if ncols > 1 else 1.0
    dy = (ydata[1]-ydata[0])/(nrows-1) if nrows > 1 else 1.0
    left, right = xdata[0]-dx/2, xdata[1]+dx/2
    top, bottom = ydata[0]-dy/2, ydata[1]+dy/2
    return [left, right, bottom, top]


matFile, outPdf, dpi = sys.argv[1], sys.argv[2], float(sys.argv[3])
data = scipy.io.loadmat(matFile, struct_as_record=False, squeeze_me=True)
axesData = data['axesData']
figW, figH = float(data['figW']), float(data['figH'])  # inches

if not isinstance(axesData, np.ndarray):
    axesData = np.array([axesData])

fig = plt.figure(figsize=(figW, figH))

for a in axesData:
    ax = fig.add_axes([a.left, a.bottom, a.width, a.height])

    fc = a.facecolor
    if as_str(fc) == 'none':
        ax.patch.set_visible(False) # transparent axes background (e.g. plotGridProps's overlay axes) - must not hide what's stacked beneath
    else:
        ax.set_facecolor(fc)

    children = a.children
    if not isinstance(children, np.ndarray):
        children = np.array([children])
    child_types = [as_str(c.type) for c in children]
    has_image = 'image' in child_types

    for c in children:
        ctype = as_str(c.type)
        if ctype == 'image':
            extent = img_extent(np.atleast_1d(c.xdata), np.atleast_1d(c.ydata), c.rgb.shape)
            alpha = np.atleast_2d(c.alpha).astype(float)
            if alpha.shape != c.rgb.shape[:2]:
                alpha = np.full(c.rgb.shape[:2], float(np.atleast_1d(c.alpha).flat[0]))
            if np.all(alpha == 1.0):
                # fully opaque - plot plain RGB with no alpha/SMask at all. Most
                # panels (e.g. plotRateMap) never need transparency; only skip this
                # for the few that do (e.g. plotGridProps's NaN-masked overlay),
                # since a PDF soft-mask is exactly the kind of construct that has
                # repeatedly caused MATLAB/Illustrator interop problems in this
                # pipeline (corrupt vector export, colour-management surprises)
                ax.imshow(c.rgb, interpolation='nearest', extent=extent)
            else:
                # bake alpha into an RGBA array rather than passing imshow's separate
                # alpha= parameter - older matplotlib (this env: 3.2.2) silently
                # ignores a 2D alpha array there (confirmed by isolated testing)
                rgba = np.dstack([c.rgb, alpha])
                ax.imshow(rgba, interpolation='nearest', extent=extent)
        elif ctype == 'line':
            x, y = np.atleast_1d(c.xdata), np.atleast_1d(c.ydata)
            ax.plot(x, y, color=c.color, linewidth=float(c.linewidth), linestyle=as_str(c.linestyle))
        elif ctype == 'area':
            x, y = np.atleast_1d(c.xdata), np.atleast_1d(c.ydata)
            base = float(c.basevalue) if np.isscalar(c.basevalue) or c.basevalue.size == 1 else 0.0
            ax.fill_between(x, base, y, facecolor=c.facecolor, edgecolor=c.edgecolor)
        elif ctype == 'text':
            ax.text(float(c.x), float(c.y), as_str(c.string), color=c.color, fontsize=float(c.fontsize),
                    ha=as_str(c.ha), va=mpl_va(as_str(c.va)), transform=ax.transData)
        elif ctype == 'figtext':
            fig.text(float(c.x), float(c.y), as_str(c.string), color=c.color, fontsize=float(c.fontsize),
                      ha=as_str(c.ha), va=mpl_va(as_str(c.va)))

    # apply MATLAB's own final axis limits/direction directly (rather than
    # guessing matplotlib's autoscale) - this is what correctly keeps e.g.
    # plotGridProps's annotation text (placed outside the image extent) in
    # view, and img_extent() above already puts the image itself in the same
    # real MATLAB data-coordinate space, so this composes correctly
    ax.set_xlim(a.xlim)
    ax.set_ylim(a.ylim)
    if as_str(a.xdir) == 'reverse':
        ax.invert_xaxis()
    if as_str(a.ydir) == 'reverse':
        ax.invert_yaxis()
    if has_image:
        ax.set_xticks([]); ax.set_yticks([])
    else:
        xtick = np.atleast_1d(a.xtick).astype(float) if np.size(a.xtick) else np.array([])
        ytick = np.atleast_1d(a.ytick).astype(float) if np.size(a.ytick) else np.array([])
        ax.set_xticks(xtick) # always set explicitly (even empty) - MATLAB may have deliberately cleared ticks
        ax.set_yticks(ytick)
        if xtick.size:
            xtl = as_str_list(a.xticklabel)
            if xtl:
                ax.set_xticklabels(xtl)
        if ytick.size:
            ytl = as_str_list(a.yticklabel)
            if ytl:
                ax.set_yticklabels(ytl)

    if as_str(a.aspectSquare):
        ax.set_aspect('equal', adjustable='box')

    xl = as_str(a.xlabel)
    yl = as_str(a.ylabel)
    tl = as_str(a.title)
    if xl:
        ax.set_xlabel(xl, fontsize=8)
    if yl:
        ax.set_ylabel(yl, fontsize=8)
    if tl:
        ax.set_title(tl, fontsize=8)
    if as_str(a.visible) == 'off':
        ax.set_frame_on(False)
        if has_image:
            ax.set_xticks([]); ax.set_yticks([])

plt.savefig(outPdf, format='pdf', dpi=dpi)
print('saved', outPdf)
