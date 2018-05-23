import dionysus as d

def test_distance():
    dgm1 = d.Diagram([(1,2), (3,4), (1., float('inf'))])
    dgm2 = d.Diagram([(0,2), (3,5), (2., float('inf'))])
    dist = d.bottleneck_distance(dgm1,dgm2)
    assert dist == 1
