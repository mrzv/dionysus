import numpy as np
import dionysus as d
import diode

points = np.random.random((100,3))

f = diode.fill_alpha_shapes(points)
f = d.Filtration(f)

r = d.homology_persistence(f, method='column_no_negative')
dgms = d.init_diagrams(r,f)
