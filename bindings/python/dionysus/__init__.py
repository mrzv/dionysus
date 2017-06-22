from ._dionysus import *

def closure(simplices, k):
    """Compute the `k`-skeleton of the closure of the list of `simplices`."""

    res = set()

    from    itertools   import combinations
    for s in simplices:
        for kk in range(1, k+2):
            for face in combinations(list(s), min(s.dimension() + 1, kk)):
                ss = Simplex(face)
                res.add(ss)

    return list(res)


def smooth(f, z, prime, show = False):
    """Smooth a given integer cocycle into a harmonic cocycle."""

    try:
        from scipy.sparse.linalg import lsqr
        from scipy.sparse        import csc_matrix
        import numpy as np
    except ImportError:
        raise ImportError("Unable to import lsqr from scipy.sparse.linalg. Have you installed scipy?")

    # Cocycle can be larger than D; we implicitly project it down

    data = []
    row  = []
    col  = []
    for i,s in enumerate(f):
        if s.dimension() == 1:
            for isb,sb in enumerate(s.boundary()):
                data.append(1. if isb % 2 == 0 else -1.)
                row.append(i)
                col.append(f.index(sb))

    z_data = [x.element if x.element < prime/2 else x.element - prime for x in z]
    z_row  = [x.index for x in z]
    z_col  = [0 for x in z]

    dim = max(max(row),max(col),max(z_row)) + 1     # max(z_row) implicitly projects the cocycle down to f
    D = csc_matrix((np.array(data), (np.array(row), np.array(col))), shape=(dim, dim))
    z = csc_matrix((z_data, (z_row, z_col)), shape=(dim, 1)).toarray()

    tol = 1e-10
    solution = lsqr(D, z, atol = tol, btol = tol, show = show)

    vertex_values = { f[i][0] : x for i,x in enumerate(solution[0]) if x != 0}
    return vertex_values
