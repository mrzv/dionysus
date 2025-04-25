import dionysus as d

w = -1      # cone vertex

def fast_zigzag(simplices, times):
    inf = float('inf')
    combined = d.LinkedMultiFiltration()
    combined.append(d.Simplex([w]), 0)
    for s,times in zip(simplices,times):
        print(s,times)
        for i,t in enumerate(times):
            if i % 2 == 0:
                sx = d.Simplex(s,t)
                combined.append(sx, len(combined))
            else:
                sx = d.Simplex(s,t).join(w)
                combined.append(sx, len(combined)-1)     # link to the previous appearance
        # if a simplex doesn't get removed, remove it at infinity
        if i % 2 == 0:
            sx = d.Simplex(s,inf).join(w)
            combined.append(sx, len(combined)-1)     # link to the previous appearance

    def cone_key(s):
        cone = w in s
        ww = cone and s.dimension() == 0
        return (not ww, cone, s.data if not cone else -s.data, s)

    combined.sort(key = cone_key)

    for i,s in enumerate(combined):
        print(combined.index(s,i), combined.index(s), s)

    return combined

def init_zigzag_diagrams(r,f):
    dgms = [d.Diagram() for _ in range(max(s.dimension() for s in f if w not in s) + 1)]

    assert w in f[0] and f[0].dimension() == 0
    for i in range(1,len(r)):
        j = r.pair(i)
        if j < i: continue      # skip negative

        if j == r.unpaired:
            j_cone = True
        else:
            j_cone = w in f[j]
        i_cone = w in f[i]

        i_data = f[i].data
        j_data = float('inf') if j == r.unpaired else f[j].data

        if not i_cone and not j_cone:
            # ordinary (closed-open)
            dgms[f[i].dimension()].append(i_data, j_data, i)
        elif i_cone and j_cone:
            # relative (open-closed)
            dgms[f[i].dimension() - 1].append(j_data, i_data, i)
        else:
            assert not i_cone and j_cone
            if i_data > j_data:     # TODO: can we check this non-numerically?
                # extended (open-open)
                dgms[f[i].dimension() - 1].append(j_data, i_data, i)
            else:
                # extended (closed-closed)
                dgms[f[i].dimension()].append(i_data, j_data, i)

    return dgms
