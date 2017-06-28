import matplotlib.pyplot as plt
from   matplotlib.colors import LogNorm, Normalize

def plot_diagram(dgm, show = False):
    """Plot the persistence diagram."""

    min_birth = min(p.birth for p in dgm)
    #max_birth = max(p.birth for p in dgm)
    #min_death = min(p.death for p in dgm)
    max_death = max(p.death for p in dgm)

    plt.axes().set_aspect('equal', 'datalim')

    plt.scatter([p.birth for p in dgm], [p.death for p in dgm])
    plt.plot([min_birth, max_death], [min_birth, max_death])        # diagonal

    ## clip the view
    #plt.axes().set_xlim([min_birth, max_birth])
    #plt.axes().set_ylim([min_death, max_death])

    if show:
        plt.show()

def plot_bars(dgm, order = 'birth', show = False):
    """Plot the barcode."""

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

    if lognorm:
        norm = LogNorm()
    else:
        norm = Normalize()

    min_birth = min(p.birth for p in dgm)
    #max_birth = max(p.birth for p in dgm)
    #min_death = min(p.death for p in dgm)
    max_death = max(p.death for p in dgm)

    plt.hist2d([p.birth for p in dgm], [p.death for p in dgm], bins = bins, norm = norm)
    plt.axes().set_aspect('equal', 'datalim')

    if diagonal:
        plt.plot([min_birth, max_death], [min_birth, max_death])        # diagonal

    ## clip the view
    #plt.axes().set_xlim([min_birth, max_birth])
    #plt.axes().set_ylim([min_death, max_death])

    plt.colorbar()

    if show:
        plt.show()
