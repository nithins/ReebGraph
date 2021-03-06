#!/usr/bin/env python
 
import sys
import vtk
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import pyrg
import attrdict


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
    
    return attrdict.AttrDict(arr=arr,shape=arr.shape,spacing=imageData.GetSpacing()) 



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