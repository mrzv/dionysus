import dionysus

# init the complex/filtration
# any filtration that produces a diagram of only infinite points will work
simplicies = [
    ([0], 0),
    ([1], 1),
]

# create the filtration
filtr = dionysus.Filtration()
for verts, idx in simplicies:
    simplex = dionysus.Simplex(verts, idx)
    filtr.append(simplex)
filtr.sort()

# create the diagram
m = dionysus.homology_persistence(filtr)
dgm = dionysus.init_diagrams(m, filtr)
dgm0 = dgm[0]

# trigger error
# line fails with error:
# Traceback (most recent call last):
#   File "/home/dave/tmp/dionysus-bug/bug.py", line 25, in <module>
#       dionysus.plot.plot_diagram(dgm0)
#   File "/.../dionysus/plot.py",
#       line 30, in plot_diagram
#       min_death = min(p.death for p in dgm if p.death != inf)
#       ValueError: min() arg is an empty sequence
dionysus.plot.plot_diagram(dgm0)

