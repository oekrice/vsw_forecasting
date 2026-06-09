"""
    Generic Python plotting routines.

    ary -- 2019 July
"""
import math
from random import shuffle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
def make_axes(
    types,
    ncols=1,
    aspect=[2, 3, 1],
    height=8,
    sharex="none",
    sharey="none",
    tex=True,
    dark=False,
):
    """
    Generate set of plot axes.
        type           -- list of plot types (1d/2d/3d) in order required for plotting
        ncols          -- number of columns in plots
        aspect         -- aspect ratios (width/height) for each type of plot
        height         -- height of figure window
        sharex, sharey -- whether to share x or y axes between subplots
    """

    # Get number of subplots from input list:
    n = len(types)
    nrows = int(math.ceil(n / float(ncols)))

    # Work out required column widths (based on top row):
    widths = [aspect[i - 1] for i in types[:ncols]]

    # Work out aspect ratio of window (width/height):
    winaspect = sum(widths[:ncols]) / float(nrows)

    # Set text to latex:
    if tex:
        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")

    # Use dark background:
    if dark:
        plt.style.use("dark_background")

    # Generate axes:
    f, ax = plt.subplots(
        nrows,
        ncols,
        gridspec_kw={"width_ratios": widths},
        figsize=(height * winaspect, height),
        sharex=sharex,
        sharey=sharey,
    )

    try:
        return ax.flatten()
    except:
        return [ax]


# ------------------------------------------------------------------------------
def list_colors(n, random=False):
    """
        Get list of RGB triplets for n colors from the Tableau20 list. https://public.tableau.com/
        Optional argument:
            random  -- whether to randomize order of resulting list
    """
    tableau20 = [
        (31, 119, 180),
        (174, 199, 232),
        (255, 127, 14),
        (255, 187, 120),
        (44, 160, 44),
        (152, 223, 138),
        (214, 39, 40),
        (255, 152, 150),
        (148, 103, 189),
        (197, 176, 213),
        (140, 86, 75),
        (196, 156, 148),
        (227, 119, 194),
        (247, 182, 210),
        (127, 127, 127),
        (199, 199, 199),
        (188, 189, 34),
        (219, 219, 141),
        (23, 190, 207),
        (158, 218, 229),
    ]

    # Scale to [0, 1] range:
    colors = []
    for i in range(n):
        if n < 10:
            # Use Tableau 10 if fewer colors, to get more contrast:
            r, g, b = tableau20[2 * i]
        else:
            r, g, b = tableau20[i % 20]
        colors.append((r / 255.0, g / 255.0, b / 255.0))

    if random:
        shuffle(colors)

    return colors


# ------------------------------------------------------------------------------
def graph(
    axes,
    x,
    f,
    xlabel="",
    ylabel="",
    title="",
    titleloc=[0.05, 0.85],
    upto=[],
    colors="",
    xlim=[],
    ylim=[],
    labels="",
    styles="",
    xaxis=False,
    yaxis=False,
    cmaplim=[0, 1],
):
    """
        Plot the graph of f against x. Both should be either numpy arrays or a list of numpy arrays.
        Optional arguments:
            xlabel, ylabel -- axis labels (strings)
            title          -- title (string)
            titleloc       -- location of the title (fraction of plot) [x, y]
            upto           -- maximum x-value to show in full weight
            colors         -- color of the plot (single string or list of strings). If it is the
                              name of a colormap, this is used to choose colors
            styles         -- linestyle of the plot (single string or list of strings)
            xlim, ylim     -- extent of plotting window [min, max] (default: full range of data]
            cmaplim        -- limits of colormap (in [0,1]) when using colormap
    """

    # If we have only one curve, put it in a list:
    if type(x) != list:
        x, f = [x], [f]
    if type(colors) != list:
        colors = [colors]
    if type(styles) != list:
        styles = [styles]
    if type(labels) != list:
        labels = [labels]
    if colors == [""]:
        colors = ["k" for x0 in x]
    if colors[0] in plt.colormaps():  # check whether this is a colormap name
        cmap = mpl.cm.get_cmap(colors[0])
        colors = cmap(np.linspace(cmaplim[0], cmaplim[1], len(x)))
    if styles == [""]:
        styles = ["-" for x0 in x]
    if styles[0] == "scatter":
        styles = ["None" for x0 in x]
        markers = ["o" for x0 in x]
    else:
        markers = ["None" for x0 in x]
    if labels == [""]:
        labels = ["" for x0 in x]
    for j in range(len(x)):
        if upto != []:
            # Plot full graph grayed out:
            axes.plot(x[j], f[j], alpha=0.2, color=colors[j], linewidth=1)
            # Plot graph up to x=upto in full weight:
            axes.plot(
                x[j][x[j] <= upto],
                f[j][x[j] <= upto],
                color=colors[j],
                linewidth=1,
                label=labels[j],
                linestyle=styles[j],
                marker=markers[j],
                markersize=2,
                markeredgewidth=0.5,
            )
        else:
            axes.plot(
                x[j],
                f[j],
                color=colors[j],
                linewidth=1,
                label=labels[j],
                linestyle=styles[j],
                marker=markers[j],
                markersize=2,
                markeredgewidth=0.5,
            )

    if xlim != []:
        xmin, xmax = xlim[0], xlim[1]
    else:
        xmin = min([x0[0] for x0 in x])
        xmax = max([x0[-1] for x0 in x])
    axes.set_xlim(xmin, xmax)

    if xaxis:
        axes.plot([xmin, xmax], [0, 0], "k-", linewidth=0.75)

    if ylim != []:
        fmin, fmax = ylim[0], ylim[1]
    else:
        fmin = min([np.min(f0) for f0 in f])
        fmax = max([np.max(f0) for f0 in f])
    axes.set_ylim(fmin, fmax)

    if yaxis:
        axes.plot([0, 0], [fmin, fmax], "k-", linewidth=0.75)

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title, x=titleloc[0], y=titleloc[1], fontsize=10)


# ------------------------------------------------------------------------------
def mesh(
    axes,
    x,
    y,
    f,
    xlabel="",
    ylabel="",
    title="",
    titleloc=[0.5, 0.85],
    upto=[],
    xlim=[],
    ylim=[],
    cmap="bwr",
    colorbar=True,
    fmax=[],
    fmin=[],
    equal=False,
    dark=False,
    noaxes=False,
    getcb=False,
    getmesh=False,
    raster=False,
):
    """
        Plot a color contour plot of f (2D numpy array) against x and y (1D numpy arrays).
    """
    # Flip dimensions if arrays are not in "meshgrid" order:
    if (np.size(f, 0) == np.size(x)) | (np.size(f, 0) == np.size(x) - 1):
        f1 = f.T
    else:
        f1 = f

    # Make background color dark:
    if dark:
        plt.gcf().patch.set_facecolor("k")

    # Colorbar limits:
    if type(fmax) == list:
        if fmax == []:
            fmax = np.max(np.abs(f))
    if type(fmin) == list:
        if fmin == []:
            fmin = -fmax

    pm = axes.pcolormesh(x, y, f1, cmap=cmap, rasterized=raster, shading="auto")
    pm.set_clim(vmin=fmin, vmax=fmax)

    # Draw line at x=upto:
    if upto != []:
        axes.plot([upto, upto], [y[0], y[-1]], "k--", linewidth=1)

    if colorbar:
        cb = plt.colorbar(pm, ax=axes)

    if xlim != []:
        axes.set_xlim(xlim[0], xlim[1])
    else:
        axes.set_xlim(x[0], x[-1])

    if ylim != []:
        axes.set_ylim(ylim[0], ylim[1])
    else:
        axes.set_ylim(y[0], y[-1])

    if equal:
        axes.set_aspect("equal")

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if dark:
        axes.set_title(title, x=titleloc[0], y=titleloc[1], color="w")
    else:
        axes.set_title(title, x=titleloc[0], y=titleloc[1])

    if noaxes:
        axes.axis("off")

    if getcb:
        return cb

    if getmesh:
        return pm
