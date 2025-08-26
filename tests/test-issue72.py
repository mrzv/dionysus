import dionysus as tpl

def detail(i,t,d,zz,cells):
    if d:
        action = 'added'
    else:
        action = 'removed'
    print(f'{i}) time= {t}, simplex:{action}')
    for z in zz:
        print(z, ' -> ', ' + '.join("%d * (%s)" % (x.element, f[cells[x.index]]) for x in z))


f_list = [[0], [1], [0,1], [2], [1,2]]
f = tpl.Filtration(f_list)
times = [[0], [1], [1, 3, 3], [2], [2]]
zz, dgms, cells = tpl.zigzag_homology_persistence(f, times, callback = detail)


f_list = [[0], [1], [0,1], [2], [1,2]]
times = [[0], [1], [1, 3, 3], [2], [2]]
cone = tpl.fast_zigzag(f_list, times)
r,v = tpl.homology_persistence(cone, method = 'matrix_v')
