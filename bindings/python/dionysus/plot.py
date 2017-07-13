def plot_diagram(dgm, show = False):
    """Plot the persistence diagram."""

    import matplotlib.pyplot as plt

    inf = float('inf')
    min_birth = min(p.birth for p in dgm if p.birth != inf)
    max_birth = max(p.birth for p in dgm if p.birth != inf)
    min_death = min(p.death for p in dgm if p.death != inf)
    max_death = max(p.death for p in dgm if p.death != inf)

    plt.axes().set_aspect('equal', 'datalim')

    min_diag = min(min_birth, min_death)
    max_diag = max(max_birth, max_death)

    plt.scatter([p.birth for p in dgm], [p.death for p in dgm])
    plt.plot([min_diag, max_diag], [min_diag, max_diag])        # diagonal

    ## clip the view
    #plt.axes().set_xlim([min_birth, max_birth])
    #plt.axes().set_ylim([min_death, max_death])

    if show:
        plt.show()

def plot_bars(dgm, order = 'birth', show = False):
    """Plot the barcode."""

    import matplotlib.pyplot as plt

    if order == 'death':
        generator = enumerate(sorted(dgm, key = lambda p: p.death))
    else:
        generator = enumerate(dgm)

    for i,p in generator:
        plt.plot([p.birth, p.death], [i,i], color = 'b')

    if show:
        plt.show()



def plot_diagram_density(dgm, bins = 200, lognorm = True, diagonal = True, show = False):
    """Plot the histogram of point density."""

    import matplotlib.pyplot as plt
    from   matplotlib.colors import LogNorm, Normalize

    if lognorm:
        norm = LogNorm()
    else:
        norm = Normalize()

    inf = float('inf')
    min_birth = min(p.birth for p in dgm if p.birth != inf)
    #max_birth = max(p.birth for p in dgm if p.birth != inf)
    #min_death = min(p.death for p in dgm if p.death != inf)
    max_death = max(p.death for p in dgm if p.death != inf)

    plt.hist2d([p.birth for p in dgm if p.birth != inf and p.death != inf], [p.death for p in dgm if p.birth != inf and p.death != inf], bins = bins, norm = norm)
    plt.axes().set_aspect('equal', 'datalim')

    if diagonal:
        plt.plot([min_birth, max_death], [min_birth, max_death])        # diagonal

    ## clip the view
    #plt.axes().set_xlim([min_birth, max_birth])
    #plt.axes().set_ylim([min_death, max_death])

    plt.colorbar()

    if show:
        plt.show()
