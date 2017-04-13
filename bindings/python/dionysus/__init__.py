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
