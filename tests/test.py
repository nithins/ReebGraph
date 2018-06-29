import pyrg
import pyrgtest as rgt
import numpy as np
import pytest

@pytest.fixture
def ds_3gauss():
    """
    make a dataset with 3 gaussian features
    """
    import numpy as np
    from scipy.stats import multivariate_normal


    def make_gaussian(size=(30,30,30),mu=(0,0,0),sigma=(0.25,0.25,0.25)):

        x, y, z = np.mgrid[-1.0:1.0:64j, -1.0:1.0:64j,-1.0:1.0:64j]
        # Need an (N, 2) array of (x, y) pairs.
        xyz = np.column_stack([x.flat, y.flat,z.flat])

        mu = np.array(mu)

        sigma = np.array(sigma)
        covariance = np.diag(sigma**2)

        f = multivariate_normal.pdf(xyz, mean=mu, cov=covariance)

        # Reshape back to a (30, 30) grid.
        f = f.reshape(x.shape)
        
        f = np.array(f,dtype=np.float32)
        
        return f

    
    arr   = make_gaussian(mu=(0,0,0))  
    arr  += make_gaussian(mu=(0.5,0,0))  
    arr  += make_gaussian(mu=(0.5,-0.5,0.5))        
    return arr

    
    
def test_disjointSets():
    rgt.testDisjointSets()


def test_toydataset():
    rgt.testToyDataset()


def test_grid3d():
    
    # 4 z slices 5in x 6 in Y
    s = (4,5,6)
    
    g  = pyrg.Grid3D(s)       
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
                
def test_compute(ds_3gauss):
    rg  = pyrg.ContourTree()
    rg.computeGrid3D(ds_3gauss)
    

def test_pickle(ds_3gauss):
    
    import cPickle as p

    rg  = pyrg.ContourTree()
    rg.computeGrid3D(ds_3gauss)
    gr = p.loads(p.dumps(rg,protocol=2))
    
    assert((rg.nodes == gr.nodes).all())
    assert((rg.arcs  == gr.arcs).all())
    assert((rg.arcmap  == gr.arcmap).all())
    
    
    
    
    
    
    
    
    
    
    
    
    
    