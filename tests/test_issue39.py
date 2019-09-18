import numpy as np
import dionysus as d

def test_issue39():
    dgm1 = np.loadtxt('data/issue39/dgm1.txt', delimiter=',')
    dgm2 = np.loadtxt('data/issue39/dgm2.txt', delimiter=',')
    dgm1 = d.Diagram(dgm1)
    dgm2 = d.Diagram(dgm2)
    dist = d.wasserstein_distance(dgm1,dgm2,q=5)
