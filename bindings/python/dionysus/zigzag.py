from collections import defaultdict
import dionysus as d

from intervaltree import IntervalTree
from sortedcontainers import SortedSet

from ._dionysus import fast_zigzag as _fast_zigzag

# This is sort of a hack; it would be better to cast whatever we get into a filtration directly in the C++ code
def fast_zigzag(f, times):
    """Build the cone to compute extended persistence equivalent to the given zigzag."""
    if type(f) is not d.Filtration:
        f = d.Filtration(f)
    return _fast_zigzag(f,times)

w = -1      # cone vertex

class ApexRepresentative:
    """Apex representative"""

    def __init__(self, dir):
        self.dir = dir
        self.horizontal = IntervalTree()
        self.vertical = []
        self.changes = SortedSet()

    def add(self, times, data):
        self.horizontal[times[0]:times[1]] = data
        self.changes.add(times[0])
        self.changes.add(times[1])

    def add_vertical(self, time, data):
        self.changes.add(time)
        self.vertical.append((time, data))

    def representative(self, time):
        if time in self.changes:
            # perturb the time
            idx = self.changes.index(time)
            if idx != 0:    # perturb left, unless at the end
                time = (self.changes[idx - 1] + self.changes[idx])/2
            else:
                time = (self.changes[idx + 1] + self.changes[idx])/2

        result = []
        for (t1,t2,(idx,c)) in self.horizontal[time]:
            result.append((idx,c))
        return result

    def __iter__(self):
        # report horizontal cells
        for (t1,t2, data) in self.horizontal:
            yield (t1,t2), data

        # report vertical cells
        for (t, data) in self.vertical:
            yield t, data

def apex(pt,r,v,f):
    """Given a point `pt` in a zigzag diagram, matrices `r` and `v`, and the cone filtration `f`, return the apex representative."""

    i = pt.data
    birth = pt.birth
    death = pt.death

    # determine the type
    j = r.pair(i)
    i_cone = w in f[i]
    j_cone = w in f[j]

    i_data = f[i].data
    j_data = f[j].data

    k = r.field()

    if not i_cone and not j_cone:
        # ordinary (closed-open)
        return lift_cycle(v[j], (death, birth), None, f, k)
    elif i_cone and j_cone:
        # relative (open-closed)
        return lift_cycle(restrict_to_base(r[j],f), (birth, death), None, f, k)
    else:
        assert not i_cone and j_cone
        if i_data > j_data:     # TODO: can we check this non-numerically?
            # extended (open-open)
            return lift_cycle(r[j], (death,birth), None, f, k)
        else:
            # extended (closed-closed)
            # There is an implicit assumption in the SoCG paper (that needs to
            # be spelled out in the journal version) that we are computing
            # H(cK,w), meaning the row of w needs to be removed.  This is
            # accounted for below, when we compute the boundary, but when r[j]
            # is 0-dimensional, it may get w in it, too, so we need to remove it.
            return lift_cycle(restrict_to_base(v[j], f), (birth, death),
                              negate(restrict_to_base(r[j],f), k), f, k)

def restrict_to_base(c, f):
    return type(c)([(x.element, x.index) for x in c if w not in f[x.index]])

# TODO: migrate this functionality to Chain
def negate(c, k):
    return type(c)([(k.neg(x.element), x.index) for x in c])

def lift_cycle(z, dir, w, fltr, k):
    # print(f"{z=},{dir=},{w=}")

    start,finish = dir

    cofaces = defaultdict(list)
    for y in z:
        s = fltr[y.index]
        if ((start < finish) and (start <= s.data <= finish)) or ((start > finish) and (start >= s.data >= finish)):
            a = k.id()    # +1
            for sb in s.boundary():
                sb_idx = fltr.index(sb, y.index)
                if sb_idx != 0:
                    cofaces[sb_idx].append((s.data, k.mul(a, y.element)))
                else:
                    # skip w; it's technically not in the boundary for H(cK,w)
                    assert sb_idx.dimension() == 0 and sb[0] == w
                a = k.neg(a)    # a = -a

    result = defaultdict(list)
    lifted = defaultdict(int)
    if w is not None:
        for x in w:
            result[x.index].append((start, x.element))
            lifted[x.index] = x.element

    # TODO: record vertical cells for completeness
    for s_idx, lst in cofaces.items():
        lst.sort(reverse = (start > finish))
        for time,coeff in lst:
            new_coeff = k.add(lifted[s_idx], coeff)
            result[s_idx].append((time,new_coeff))
            lifted[s_idx] = new_coeff

    for idx,coeff in lifted.items():
        if not k.is_zero(coeff):
            result[idx].append((finish, 0))

    for i,lst in result.items():
        lst.sort()

    # build interval trees
    rep = ApexRepresentative(dir)
    for idx, lst in result.items():
        for (t1,c),(t2,_) in zip(lst, lst[1:]):
            if t1 != t2:
                rep.add((t1,t2), (idx,c))

    return rep

def point_representative(apex_representative, time):
    """Given an apex representative, return zigzag representative at the given `time`."""

    return apex_representative.representative(time)
