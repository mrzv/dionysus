import dionysus as d

simplices = [[0], [1], [0,1], [2], [0,2], [1,2]]
times = [[.4,.6,.7,1.], [.1,.2,.3,1.], [.8,.95], [.5], [.8,1.], [.9,1.]]

for s,ts in zip(simplices,times):
    print(s, ts)

cone = d.fast_zigzag(simplices, times)
r,v = d.homology_persistence(cone, method = 'matrix_v')

dgms = d.init_zigzag_diagrams(r,cone)
print(dgms)

max_t = max(max(t) for t in times)
for dim,type_dgm in enumerate(dgms):
    print("Dimension:", dim)
    for t,dgm in type_dgm.items():
        print("Type:", t)
        for pt in dgm:
            print(pt)
            apex_rep = d.apex(pt,r,v,cone)
            print("apex representative: ", end='')
            for (t1,t2, (x,c)) in apex_rep:
                print(f"{cone[x]} × [{t1},{t2}] ⋅ {c}")

            if pt.death != float('inf'):
                middle = (pt.birth + pt.death)/2
            else:
                middle = max_t + 1
            middle_representative = d.point_representative(apex_rep, middle)
            print(f"midpoint ({middle}) representative:", ' + '.join(f"{coeff} ⋅ {cone[idx]}" for (idx,coeff) in middle_representative))
