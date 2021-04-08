import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


new_rc_params = {
    # Default parameters for plot
    #* text
    # 'text.usetex': True,
    'font.size': 50,
    # "text.latex.preamble": r"\usepackage{amsmath}",

    #* figure
    'figure.figsize': (10, 10),
    'figure.facecolor': 'white',
    'axes.linewidth': 3,
    'image.aspect': 'equal',

    #* line
    'lines.linewidth': 4,
    'lines.markersize': 15,

    #* legend
    'legend.frameon': False,
    'legend.fontsize': 40,
    'legend.handlelength': 1.0,
    'legend.handletextpad': 0.5,
    'legend.labelspacing': 0.2,
    'legend.borderaxespad': 0.4,

    #* label
    'axes.labelsize': 50,

    #* tick
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 3,
    'ytick.major.width': 3,
    'xtick.major.size': 10,
    'ytick.major.size': 10,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'xtick.minor.size': 5,
    'ytick.minor.size': 5,

    #* tick label
    'xtick.labelsize': 36,
    'ytick.labelsize': 36,
    'xtick.major.pad': 15,
    'ytick.major.pad': 15,

    #* save
    'savefig.dpi': 300,
    # 'savefig.transparent': True,
    'savefig.facecolor': "#ffffff",
    'savefig.bbox': 'tight',
    'savefig.format': 'pdf'
}

mpl.rcParams.update(new_rc_params)

#* log-log scale linear fitting of (x,y). output two points at x=xFit[0], xFit[1] including offset
#* return xFit, yFit, gradient, residual


def logFit(x, y, start=None, end=None, offset=0.0):
    if not start:
        start = x[0]
    if not end:
        end = x[-1]
    poly, residual, _, _, _ = np.polyfit(np.log10(x), np.log10(y), 1, full=True)
    fitX = np.array([start, end])
    fitY = np.power(10.0, poly[1] - offset) * np.power(fitX, poly[0])
    return fitX, fitY, poly[0], residual


def linLogFit(x, y, start=None, end=None, offset=0.0):
    if not start:
        start = x[0]
    if not end:
        end = x[-1]
    poly, residual, _, _, _ = np.polyfit(x, np.log10(y), 1, full=True)
    fitX = np.array([start, end])
    fitY = np.power(10.0, poly[0] * fitX + poly[1] - offset)
    return fitX, fitY, poly[0], residual


def logLinFit(x, y, start=None, end=None, offset=0.0):
    if not start:
        start = x[0]
    if not end:
        end = x[-1]
    poly, residual, _, _, _ = np.polyfit(np.log10(x), y, 1, full=True)
    fitX = np.array([start, end])
    fitY = poly[0] * np.log10(fitX) + poly[1] - offset
    return fitX, fitY, poly[0], residual


def loglogLine(x0, y0, slope, x1, x2=None, offset=0.0):
    assert y0 != offset, "At log-log-line, check the y0 and offset"
    firstY = np.power(x1 / x0, slope) * (y0 - offset)
    if x2:
        secondY = np.power(x2 / x0, slope) * (y0 - offset)
    else:
        secondY = y0 - offset
    return firstY, secondY


def loglog_trans_y(x1, y1, x2, y2, offset=0.0):
    newX = x1, x2
    newY = y1 * np.power(10.0, -offset), y2 * np.power(10.0, -offset)
    return newX, newY


def loglog_trans_x(x1, y1, x2, y2, offset=0.0):
    newX = x1 * np.power(10.0, -offset), x2 * np.power(10.0, -offset)
    newY = y1, y2
    return newX, newY


if __name__ == "__main__":
    print("This is a module draw.py")
