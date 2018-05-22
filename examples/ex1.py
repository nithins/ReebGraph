import pyrg as rg
import numpy as np
from scipy.stats import multivariate_normal


def make_gaussian(size=(30,30,30),mu=(0,0,0),sigma=(0.25,0.25,0.25)):

    x, y, z = np.mgrid[-1.0:1.0:30j, -1.0:1.0:30j,-1.0:1.0:30j]
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


def make_function():
    #make a gaussian with 3 features
    arr   = make_gaussian(mu=(0,0,0))  
    arr  += make_gaussian(mu=(0.5,0,0))  
    arr  += make_gaussian(mu=(0.5,-0.5,0.5))
    return arr



    

if __name__=="__main__":
    #make_function()
    
    #g  = rg.Grid3D((30,30,30))
    #a  = np.array(g, copy = False)    
    #a[:] = make_function()
    
    #mt =  rg.MergeTree()
    #mt.computeContourTree(g)
    #s,w = mt.genPersistenceHierarchy()
    
    #print s
    #print w
    
    nodes,arcs,amap =  rg.computeCT_Grid3D(make_function())
    aseq,awts,fh  = rg.simplifyCT_Pers(nodes,arcs)

    
    print
    print nodes
    print arcs    
    print aseq
    print awts
    print fh
