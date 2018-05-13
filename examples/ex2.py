#!/usr/bin/env python
 
import sys
import vtk
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
 
class MainWindow(QtGui.QMainWindow):
 
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)
 
        self.frame = QtGui.QFrame()
 
        self.vl = QtGui.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)
 
        self.ren = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
         
        # Create Dataset
        self.dataset = self.create_dataset();
        
        # Create Volume
        volume = self.create_volume()
 
        self.ren.AddVolume(volume)
 
        self.ren.ResetCamera()
 
        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)
 
        self.show()
        self.iren.Initialize()
        
    def create_dataset(self):
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

        #make a gaussian with 3 features
        arr   = make_gaussian(mu=(0,0,0))  
        arr  += make_gaussian(mu=(0.5,0,0))  
        arr  += make_gaussian(mu=(0.5,-0.5,0.5))        
        return arr
        
        #data_matrix = np.zeros([75, 75, 75], dtype=np.uint8)
        #data_matrix[0:35, 0:35, 0:35] = 50
        #data_matrix[25:55, 25:55, 25:55] = 100
        #data_matrix[45:74, 45:74, 45:74] = 150

        #return data_matrix
        
    def create_volume(self):
        
        from vtk.util import numpy_support as nps
        
        # We begin by creating the data we want to render.
        # For this tutorial, we create a 3D-image containing three overlaping cubes. 
        # This data can of course easily be replaced by data from a medical CT-scan or anything else three dimensional.
        # The only limit is that the data must be reduced to unsigned 8 bit or 16 bit integers.
        
        arr = self.dataset
        fmin,fmax = np.min(arr),np.max(arr)
        
        def nfv(t):
            return fmin+t*(fmax - fmin)
                        
        
        scalars = nps.numpy_to_vtk(arr.ravel())
        scalars.SetName("Scalars")
        
        imageData = vtk.vtkImageData()
        
        imageData.SetDimensions(arr.shape)
        #assume 0,0 origin and 1,1 spacing.
        #__depthImageData.SetSpacing([1,1])
        #__depthImageData.SetOrigin([0,0])
        imageData.GetPointData().SetScalars(scalars)
        imageData.SetExtent(0, arr.shape[2]-1, 0, arr.shape[1]-1, 0, arr.shape[0]-1)
        imageData.SetWholeExtent(0, arr.shape[2]-1, 0, arr.shape[1]-1, 0, arr.shape[0]-1)
        imageData.Update()  

        # The following class is used to store transparencyv-values for later retrival. In our case, we want the value 0 to be
        # completly opaque whereas the three different cubes are given different transperancy-values to show how it works.
        alphaChannelFunc = vtk.vtkPiecewiseFunction()
        alphaChannelFunc.AddPoint(nfv(0.0) ,  0.0)
        alphaChannelFunc.AddPoint(nfv(0.2) ,  0.01)
        alphaChannelFunc.AddPoint(nfv(0.5)  ,  0.1)
        alphaChannelFunc.AddPoint(nfv(1.0)  , 0.2)

        # This class stores color data and can create color tables from a few color points. For this demo, we want the three cubes
        # to be of the colors red green and blue.
        colorFunc = vtk.vtkColorTransferFunction()
        colorFunc.AddRGBPoint(nfv(0.01) , 0.0, 0.0, 1.0)
        colorFunc.AddRGBPoint(nfv(0.5)  , 1.0, 1.0, 1.0)
        colorFunc.AddRGBPoint(nfv(1.0)  , 1.0, 0.0, 0.0)

        # The preavius two classes stored properties. Because we want to apply these properties to the volume we want to render,
        # we have to store them in a class that stores volume prpoperties.
        volumeProperty = vtk.vtkVolumeProperty()
        volumeProperty.SetColor(colorFunc)
        volumeProperty.SetScalarOpacity(alphaChannelFunc)

        # This class describes how the volume is rendered (through ray tracing).
        compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
        # We can finally create our volume. We also have to specify the data for it, as well as how the data will be rendered.
        volumeMapper = vtk.vtkOpenGLGPUVolumeRayCastMapper()
        #volumeMapper.SetVolumeRayCastFunction(compositeFunction)
        volumeMapper.SetInputConnection(imageData.GetProducerPort())

        # The class vtkVolume is used to pair the preaviusly declared volume as well as the properties to be used when rendering that volume.
        volume = vtk.vtkVolume()
        volume.SetMapper(volumeMapper)
        volume.SetProperty(volumeProperty)
        
        return volume

 
 
if __name__ == "__main__":
 
    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow()
 
    sys.exit(app.exec_())