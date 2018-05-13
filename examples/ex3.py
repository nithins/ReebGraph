import os,sys

from PyQt4.QtCore import *
from PyQt4 import  QtGui
from PyQt4.QtWebKit import *

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import pyrg
import vtk


def create_dataset():
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
        
        # add a renderer to the widget
        self.ren = vtk.vtkRenderer()
        self.renderWindow.AddRenderer(self.ren)
        
        # add a volume and ResetCamera
        self.ren.AddVolume(volume) 
        self.ren.ResetCamera()
        
        #prepare interactor
        self.iren = self.renderWindow.GetInteractor()
        self.iren.Initialize()



class ReebgraphModel(QObject):
    @pyqtSlot()
    def quit(self):
        QApplication.quit()

class WebViewWindow(QtGui.QWidget):
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
#        super(Window, self).__init__()
        self.view = QWebView(self)

        self.setupInspector()

        self.splitter = QtGui.QSplitter(self)
        self.splitter.setOrientation(Qt.Vertical)

        layout = QtGui.QVBoxLayout(self)
        layout.setMargin(0)
        layout.addWidget(self.splitter)

        self.splitter.addWidget(self.view)
        self.splitter.addWidget(self.webInspector)
        
        
        shortcut = QtGui.QShortcut(self)
        shortcut.setKey(Qt.Key_F5)
        shortcut.activated.connect(self.refresh)
        
        self.rgm = ReebgraphModel(self)
        self.view.page().mainFrame().addToJavaScriptWindowObject("rgm", self.rgm)
        
        self.refresh()


        
    def refresh(self):
        self.view.setHtml(open("force.html").read())

    def setupInspector(self):
        page = self.view.page()
        page.settings().setAttribute(QWebSettings.DeveloperExtrasEnabled, True)
        self.webInspector = QWebInspector(self)
        self.webInspector.setPage(page)

        shortcut = QtGui.QShortcut(self)
        shortcut.setKey(Qt.Key_F12)
        shortcut.activated.connect(self.toggleInspector)
        self.webInspector.setVisible(False)

    def toggleInspector(self):
        self.webInspector.setVisible(not self.webInspector.isVisible())
        
class MainWindow(QtGui.QMainWindow):
 
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)
        
        # Create Dataset
        self.dataset = create_dataset();   
        
        # Compute Reeb graph
        self.rg = pyrg.computeCT_Grid3D(self.dataset)        
        
        # Create UI 
        
        self.splitter = QtGui.QSplitter(self)
        self.splitter.setOrientation(Qt.Horizontal)
        
        self.vtkWidget_vr = QVTKRenderWindowInteractor(self.splitter)
        self.webview = WebViewWindow(self.splitter)
 
        self.layout = QtGui.QHBoxLayout()
        self.layout.setMargin(0)
        self.layout.addWidget(self.splitter)
        
        self.splitter.addWidget(self.vtkWidget_vr)
        self.splitter.addWidget(self.webview)

        
        self.setCentralWidget(self.splitter)
        
        self.show()          
        
        ## Create Volume Render pipeline
        self.volumeRenderPipeline = VolumeRenderPipeine(self.dataset,self.vtkWidget_vr.GetRenderWindow()) 
 
        ## Create Reebgraph Render pipeline
        #self.reebgraphRenderPipeline = ReebgraphRenderPipeline(self.dataset,self.vtkWidget_rg.GetRenderWindow()) 


def main():
    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow()
 
    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()
