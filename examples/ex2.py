#!/usr/bin/env python
 
import sys
import vtk
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import pyrg
#import attrdict

def create_3gauss():
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

    #data_matrix = np.zeros([75, 75, 75], dtype=np.uint8)
    #data_matrix[0:35, 0:35, 0:35] = 50
    #data_matrix[25:55, 25:55, 25:55] = 100
    #data_matrix[45:74, 45:74, 45:74] = 150

    #return data_matrix
    
def read_vti(f):
    
    from vtk.util import numpy_support as nps
    
    reader = vtk.vtkXMLImageDataReader() if f.endswith(".vti") else vtk.vtkStructuredPointsReader()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    reader.SetFileName(f);
    reader.Update();   
    
    
    imageData = reader.GetOutput()    
    
    dim = tuple(reversed(imageData.GetDimensions()))    
    arr = nps.vtk_to_numpy(imageData.GetPointData().GetArray(0))
    arr = np.array(arr,np.float32).reshape(dim)
    arr = (arr - arr.min())/(arr.max() - arr.min())
    
    #return attrdict.AttrDict(arr=arr,shape=arr.shape,spacing=imageData.GetSpacing())
    return arr



class VolumeRenderPipeine:
    """
    Given an input numpy volume this will volume render to given vtkRenderWindow
    """ 
    def __init__(self, arr, renderWindow):
        self.arr          = arr          # input volume
        self.renderWindow = renderWindow #output render window
        
        
        from vtk.util import numpy_support as nps
        
        # We begin by creating the data we want to render.
        # For this tutorial, we create a 3D-image containing three overlaping cubes. 
        # This data can of course easily be replaced by data from a medical CT-scan or anything else three dimensional.
        
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
        volumeProperty.ShadeOn()
        volumeProperty.SetScalarOpacity(alphaChannelFunc)

        # We can finally create our volume. We also have to specify the data for it, as well as how the data will be rendered.
        volumeMapper = vtk.vtkOpenGLGPUVolumeRayCastMapper()        
        volumeMapper.SetInputData(imageData) if vtk.VTK_MAJOR_VERSION > 5 else volumeMapper.SetInputConnection(imageData.GetProducerPort())

        # The class vtkVolume is used to pair the preaviusly declared volume as well as the properties to be used when rendering that volume.
        volume = vtk.vtkVolume()
        volume.SetMapper(volumeMapper)
        volume.SetProperty(volumeProperty)
        
        # Add a bounding box around the dataset
        bbFilter = vtk.vtkOutlineFilter()
        bbFilter.SetInputData(imageData) if vtk.VTK_MAJOR_VERSION > 5 else bbFilter.SetInputConnection(imageData.GetProducerPort())

        bbMapper = vtk.vtkDataSetMapper()
        bbMapper.SetInputConnection(bbFilter.GetOutputPort())

        bbActor = vtk.vtkActor()
        bbActor.GetProperty().EdgeVisibilityOn()
        bbActor.GetProperty().SetEdgeColor(1,1,1)
        bbActor.SetMapper(bbMapper)

        
        # add a renderer to the widget
        self.ren = vtk.vtkRenderer()
        self.renderWindow.AddRenderer(self.ren)
        
        # add a volume and ResetCamera
        self.ren.AddVolume(volume) 
        self.ren.AddActor(bbActor)
        self.ren.ResetCamera()
        
        #prepare interactor
        istyle = vtk.vtkInteractorStyleTrackballCamera()
        self.iren = self.renderWindow.GetInteractor()
        self.iren.SetInteractorStyle(istyle)
        self.iren.Initialize()

        self.imageData = imageData
        self.volumeMapper = volumeMapper
        
        
    def reloadData(self):
        self.imageData.Modified()
        #self.imageData.Update()
        self.volumeMapper.Update()
        self.renderWindow.Render()
        
    def setSpacing(self,v):
        self.imageData.SetSpacing(v)
        self.reloadData()


class ReebgraphRenderPipeline:
    """
    Given an input numpy volume this will volume render to given vtkRenderWindow
    """ 
    def __init__(self, rg, renderWindow):
        
        self.rg           = rg            # input reeb graph
        self.renderWindow = renderWindow  # output render window
        
        # Create source
        source = vtk.vtkSphereSource()
        source.SetCenter(0, 0, 0)
        source.SetRadius(5.0)
 
        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
 
        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
         
        
        # add a renderer to the widget
        self.ren = vtk.vtkRenderer()
        self.renderWindow.AddRenderer(self.ren)
        
        # add a actors and ResetCamera
        self.ren.AddActor(actor) 
        self.ren.ResetCamera()
        
        #prepare interactor
        self.iren = self.renderWindow.GetInteractor()
        self.iren.Initialize()
        

 
class MainWindow(QtGui.QMainWindow):
 
    def __init__(self, parent = None, dataset=None):
        QtGui.QMainWindow.__init__(self, parent)
        
        # Create Dataset
        self.dataset = create_3gauss() if dataset is None else dataset
        
        # Compute Reeb graph
        #self.rg = pyrg.computeCT_Grid3D(self.dataset)        
        
        # Create UI 
        self.frame = QtGui.QFrame()
 
        self.hl = QtGui.QHBoxLayout()
        self.vtkWidget_vr = QVTKRenderWindowInteractor(self.frame)
        self.vtkWidget_rg = QVTKRenderWindowInteractor(self.frame)
        self.hl.addWidget(self.vtkWidget_vr)
        self.hl.addWidget(self.vtkWidget_rg)
        
        self.frame.setLayout(self.hl)
        self.setCentralWidget(self.frame)
 
        self.show()          
        
        # Create Volume Render pipeline
        self.volumeRenderPipeline = VolumeRenderPipeine(self.dataset,self.vtkWidget_vr.GetRenderWindow()) 
 
        # Create Reebgraph Render pipeline
        self.reebgraphRenderPipeline = ReebgraphRenderPipeline(self.dataset,self.vtkWidget_rg.GetRenderWindow()) 
        
        
                


 
 
if __name__ == "__main__":
 
    
    dataset = None
    
    if len(sys.argv) >1:
        dataset = read_vti(sys.argv[1])

    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow(dataset=dataset)
 
    sys.exit(app.exec_())
        