import dionysus as d
import matplotlib.pyplot as plt

def test_issue47():
    # init the complex/filtration
    # any filtration that produces a diagram of only infinite points will work
    simplicies = [
        ([0], 0),
        ([1], 1),
    ]

    # create the filtration
    filtr = d.Filtration()
    for verts, idx in simplicies:
        simplex = d.Simplex(verts, idx)
        filtr.append(simplex)
    filtr.sort()

    # create the diagram
    m = d.homology_persistence(filtr)
    dgm = d.init_diagrams(m, filtr)

    d.plot.plot_diagram(dgm[0])
