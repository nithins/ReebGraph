import pyrg as rg
import pyrgtest as rgt
import numpy as np
    
    
def test_disjointSets():
    rgt.testDisjointSets()


def test_toydataset():
    rgt.testToyDataset()


def test_grid3d():
    
    # 4 z slices 5in x 6 in Y
    s = (4,5,6)
    
    g  = rg.Grid3D(s)       
    a  = np.array(g, copy = False)
    a[:] = np.arange(s[0]*s[1]*s[2]).reshape(s)
    
    assert a.shape == s
    print s 
    print a
    print g
    
    for z in range(s[0]):
        for y in range(s[1]):
            for x in range(s[2]):
                assert a[z,y,x] == g.atXYZ(x,y,z)     
    

    
