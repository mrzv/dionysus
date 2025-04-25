import dionysus as d

simplices = [[0], [1], [0,1], [2], [0,2], [1,2]]
times = [[.4,.6,.7,1.], [.1,.2,.3,1.], [.8,.95], [.5], [.8,1.], [.9,1.]]

cone = d.fast_zigzag(simplices, times)
r,v = d.homology_persistence(cone, method = 'matrix_v')

dgms = d.init_zigzag_diagrams(r,cone)
print(dgms)

# TODO: lift apex representatives

for dim,dgm in enumerate(dgms):
    print("Dimension:", dim)
    for pt in dgm:
        print(pt)
