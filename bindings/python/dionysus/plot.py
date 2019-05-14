def plot_diagram(dgm, show=False, labels=False, ax=None,
                 line_style=None, pt_style=None):
    """
    Plot the persistence diagram.

    Arguments:
        dgm (Diagram): See for example `init_diagrams`.

    Keyword Arguments:
        show (bool): Display the plot. (Default: False)
        labels (bool): Set axis labels. (Default: False)
        ax (AxesSubplot): Axes that should be used for plotting (Default: None)
        pt_style (dict): argments passed to `ax.scatter` for style of points.
        line_style (dict): argments passed to `ax.plot` for style of diagonal line.
    """

    import matplotlib.pyplot as plt

    line_kwargs = {}
    pt_kwargs = {}
    if pt_style is not None:
        pt_kwargs.update(pt_style)
    if line_style is not None:
        line_kwargs.update(line_style)


    inf = float('inf')
    min_birth = min(p.birth for p in dgm if p.birth != inf)
    max_birth = max(p.birth for p in dgm if p.birth != inf)
    min_death = min(p.death for p in dgm if p.death != inf)
    max_death = max(p.death for p in dgm if p.death != inf)

    if ax is None:
        ax = plt.axes()
    ax.set_aspect('equal', 'datalim')

    min_diag = min(min_birth, min_death)
    max_diag = max(max_birth, max_death)
    ax.scatter([p.birth for p in dgm], [p.death for p in dgm], **pt_kwargs)
    ax.plot([min_diag, max_diag], [min_diag, max_diag], **line_kwargs)

    if labels:
        ax.set_xlabel('birth')
        ax.set_ylabel('death')

    ## clip the view
    #plt.axes().set_xlim([min_birth, max_birth])
    #plt.axes().set_ylim([min_death, max_death])

    if show:
        plt.show()

def plot_bars(dgm, order='birth', show=False, ax=None, **bar_style):
    """
    Plot the barcode.

    Arguments:
        dgm (Diagram): See for example `init_diagrams`.

    Keyword Arguments:
        order (str): How to sort the bars, either 'death' or 'birth'
                     (Default: 'birth')
        show (bool): Display the plot. (Default: False)
        ax (AxesSubplot): Axes that should be used for plotting (Default: None)
        **bar_style: Arguments passed to `ax.plot` for style of the bars.
                     (Defaults: color='b')
    """

    import matplotlib.pyplot as plt

    bar_kwargs = {'color': 'b'}
    bar_kwargs.update(bar_style)

    if order == 'death':
        generator = enumerate(sorted(dgm, key = lambda p: p.death))
    else:
        generator = enumerate(dgm)

    if ax is None:
        ax = plt.axes()

    for i,p in generator:
        ax.plot([p.birth, p.death], [i,i], **bar_kwargs)

    if show:
        plt.show()


def plot_diagram_density(dgm, lognorm=True, diagonal=True,
                         show=False, labels=False, ax=None, **hist_style):
    """
    Plot the histogram of point density.

    Arguments:
        dgm (Diagram): See for example `init_diagrams`.

    Keyword Arguments:
        bins (int): bins for histogram, see `ax.hist2d` (Default: 200)
        lognorm (bool): Use logarithmic norm (Default: True)
        diagonal (bool):  (Default: True)
        show (bool): Display the plot. (Default: False)
        labels (bool): Set axis labels. (Default: False)
        ax (AxesSubplot): Axes that should be used for plotting (Default: None)
        **hist_style: Arguments passed to `ax.hist2d` for style of the histogram.
            (Defaults: bins=200)
    """

    import matplotlib.pyplot as plt
    from   matplotlib.colors import LogNorm, Normalize

    hist_kwargs = {'bins': 200}
    hist_kwargs.update(hist_style)

    if lognorm:
        norm = LogNorm()
    else:
        norm = Normalize()

    inf = float('inf')
    min_birth = min(p.birth for p in dgm if p.birth != inf)
    #max_birth = max(p.birth for p in dgm if p.birth != inf)
    #min_death = min(p.death for p in dgm if p.death != inf)
    max_death = max(p.death for p in dgm if p.death != inf)

    if ax is None:
        _, ax = plt.subplots()

    hist2d, histx, histy, im = ax.hist2d([p.birth for p in dgm if p.birth != inf and p.death != inf], [p.death for p in dgm if p.birth != inf and p.death != inf], norm=norm, **hist_kwargs)
    ax.set_aspect('equal', 'datalim')
    if labels:
        ax.set_xlabel('birth')
        ax.set_ylabel('death')

    if diagonal:
        ax.plot([min_birth, max_death], [min_birth, max_death])

    ## clip the view
    #plt.axes().set_xlim([min_birth, max_birth])
    #plt.axes().set_ylim([min_death, max_death])

    if labels:
        plt.colorbar(im, ax=ax, label='overlap quantity')
    else:
        plt.colorbar(im, ax=ax)

    if show:
        plt.show()
