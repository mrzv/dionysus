import numpy as np
import dionysus as d
import diode

points = np.random.random((100,3))

f = diode.fill_alpha_shapes(points)
f = d.Filtration(f)

r,v = d.homology_persistence(f, method="matrix_v")
dgms = d.init_diagrams(r,f)

for i,(x,y) in enumerate(zip(r,v)):
    print(f"R[{i}] =", x)
    print(f"V[{i}] =", y)
