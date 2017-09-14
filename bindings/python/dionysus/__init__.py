from __future__ import absolute_import
from ._dionysus import *
from . import plot     # make it available without an explicit import
from ._version import version_info, __version__

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

    row_max = max(row)           # x.index <= row_max condition below projects the cocycle to the filtration
    z_data = [x.element if x.element < prime/2 else x.element - prime for x in z if x.index <= row_max]
    z_row  = [x.index for x in z if x.index <= row_max]
    z_col  = [0 for x in z if x.index <= row_max]

    dim = max(row_max,max(col)) + 1
    D = csc_matrix((np.array(data), (np.array(row), np.array(col))), shape=(dim, dim))
    z = csc_matrix((z_data, (z_row, z_col)), shape=(dim, 1)).toarray()

    tol = 1e-10
    solution = lsqr(D, z, atol = tol, btol = tol, show = show)

    max_vrt = max(s[0] for s in f if s.dimension() == 0)

    vertex_values = [0. for _ in range(max_vrt + 1)]
    for i,x in enumerate(solution[0]):
        if f[i].dimension() == 0:
            vertex_values[f[i][0]] = x

    return vertex_values
