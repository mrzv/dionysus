import dionysus as d

simplices = [[0], [1], [0,1], [2], [0,2], [1,2]]
times = [[.4,.6,.7,1.], [.1,.2,.3,1.], [.8,.95], [.5], [.8,1.], [.9,1.]]

cone = d.fast_zigzag(simplices, times)
r,v = d.homology_persistence(cone, method = 'matrix_v')

dgms = d.init_zigzag_diagrams(r,cone)
print(dgms)

# TODO: lift apex representatives

for dim,type_dgm in enumerate(dgms):
    print("Dimension:", dim)
    for t,dgm in type_dgm.items():
        print("Type:", t)
        for pt in dgm:
            print(pt)
            apex_rep = d.apex(pt,r,v,cone)
            for x,lst in apex_rep.items():
                print(f"{cone[x]}: ", end='')
                for (t1,c1),(t2,c2) in zip(lst,lst[1:]):
                    print(f"[{t1},{t2}] ⋅ {c1}", end = '  ')
                print()
