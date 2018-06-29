import os,sys

from PyQt4.QtCore import *
from PyQt4 import  QtGui
from PyQt4.QtWebKit import *

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import pyrg
import vtk

from ex2 import create_3gauss,VolumeRenderPipeine,read_vti
import json,attrdict
import cPickle as pickle


class ReebgraphModel(QObject):
    """
    Class that handles creating and maintaining the reebgraph
    Also provides jsonified data for webview
    """
    
    dsChanged = pyqtSignal()    

    def __init__(self, dataset,rg, parent = None):
        QObject.__init__(self, parent)
        self.ods = np.copy(dataset)
        self.ds  = dataset        
        self.rg  = rg
        
        self.arcTree = self.make_arcTree()        

    def get_type(self,a,b):
        """
        Returns the type of arc a,b. min-sad, sad-sad, sad-max. 
        """
        return ("min" if self.rg.nodes[a]["type"] == 1 else "sad") +"-" + ("max" if self.rg.nodes[b]["type"] == 2 else "sad")

    def make_arcTree(self) :
        
        arcTree = {}
        
        for a,b in self.rg.arcs:
            arcTree[a,b] = {
                "name":str((a,b)),
                "pers":float(self.rg.nodes[b]["fn"] - self.rg.nodes[a]["fn"]),
                "type":self.get_type(a,b),
                "weight":1.0,
                "par": None,
                "selected":False,
                }
                                 
        for (t,c,m,a,b),w in zip(self.rg.fhier,self.rg.fwts):
            
            arcTree[(a,b)] = {
                "name":str((a,b)),
                "pers":float(self.rg.nodes[b]["fn"] - self.rg.nodes[a]["fn"]),
                "type":self.get_type(a,b),
                "weight":1.0,
                "par": None,
                "selected":False,
                }

            clu = (c,m) if t == 0 else (m,c)
            
            arcTree[clu]["par"] = (a,b); arcTree[clu]["weight"] = float(w)
            arcTree[m,b]["par"] = (a,b); arcTree[m,b]["weight"] = float(w)
            arcTree[a,m]["par"] = (a,b); arcTree[a,m]["weight"] = float(w)
        
        return arcTree

    def adjust_ds(self):
                    
        arcRemap = np.zeros(len(self.rg.arcs),np.float32)
        
        for ano in range(len(self.rg.arcs)):
            a,b = tuple(map(int,self.rg.arcs[ano]))
            
            n  = self.arcTree[a,b]
            s  = n["selected"]
            while n["par"] != None:
                n  = self.arcTree[n["par"]]
                s |= n["selected"]
            arcRemap[ano] = s
            
        self.ds[:] = self.ods*arcRemap[self.rg.arcmap]
        self.dsChanged.emit()


    
    @pyqtSlot(result=str)
    def json(self):
        
        nodes,arcs = self.rg.nodes,self.rg.arcs

        rng   = [float(nodes["fn"].min()),float(nodes["fn"].max())]                
        nodes = [ {"id":i, "name":str(n["id"]), "fn":float(n["fn"]),"group":int(n["type"]), "weight":1.0 } for i,n in enumerate(nodes)]
        links = [ {"source":arc[0], "target":arc[1]} for arc in arcs ]
        
        for (a,b),w in zip(self.rg.farcs,self.rg.fwts):
            nodes[a]["weight"] = float(w);
            nodes[b]["weight"] = float(w);

        
        return json.dumps({"nodes":nodes,"links":links,"range":rng})

    @pyqtSlot(result=str)
    def hierTree_json(self):
        
        nodes,arcs = self.rg.nodes,self.rg.arcs
        
        arcTreeTD = {}
        
        for a,b in self.rg.arcs:
            a,b = int(a),int(b)
                        
            arcTreeTD[str((a,b))] = {
                "name":str((a,b)),
                #"size":float(nodes[b]["fn"] - nodes[a]["fn"]),
                #"size":10,
                "size":1 + 100*float(self.rg.nodes[b]["fn"] - self.rg.nodes[a]["fn"]),
                "type":self.get_type(a,b),
                "weight":1.0,
                "selected":self.arcTree[(a,b)]["selected"],
                }
        
        
                                 
        for (t,c,m,d,u),w in zip(self.rg.fhier,self.rg.fwts):
                        
            clu = str((c,m) if t == 0 else (m,c)); m_u = str((m,u));  d_m = str((d,m)); d_u = str((d,u))
            
            arcTreeTD[clu]["weight"] = arcTreeTD[m_u]["weight"] = arcTreeTD[d_m]["weight"] = float(w)
            
            assert not arcTreeTD.has_key((d_u))
            arcTreeTD[d_u] = {
                "name":d_u,
                "children": [arcTreeTD[clu],arcTreeTD[d_m],arcTreeTD[m_u]],
                #"size":float(nodes[u]["fn"] - nodes[d]["fn"]),
                "type":self.get_type(d,u),
                "weight":1.0,                
                "selected":self.arcTree[int(d),int(u)]["selected"],
                }
            
            del arcTreeTD[m_u]; del arcTreeTD[d_m]; del arcTreeTD[clu]
                                
        for k,v in arcTreeTD.iteritems():
            arcTreeTD = v
            break
        
        return json.dumps(arcTreeTD)

    @pyqtSlot(str,bool,result=bool)
    def selectArc(self,aname,selected):
        try:
            a,b = tuple(map(int,aname[1:-1].split(",")))
            self.arcTree[a,b]["selected"] = selected
            self.adjust_ds()
        except Exception as e:
            print "Selection Failed aname=",aname,"e=",e
            return False
        return True
            
        

    
    @pyqtSlot()
    def quit(self):
        QApplication.quit()
        

class WebViewWindow(QtGui.QWidget):
    def __init__(self,parent = None):
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
        
        self.rgm = None
        
        
        shortcut = QtGui.QShortcut(self)
        shortcut.setKey(Qt.Key_F5)
        shortcut.activated.connect(self.refresh)
        
        self.view.page().mainFrame().javaScriptWindowObjectCleared.connect(self.addJsObjects)
                

    def addJsObjects(self):
        if self.rgm != None:
            self.view.page().mainFrame().addToJavaScriptWindowObject("rgm", self.rgm)
        if self.app != None:
            self.view.page().mainFrame().addToJavaScriptWindowObject("app", self.app)
        
        
    def refresh(self,rgm=None,app=None):
        
        if rgm is not None:
            self.rgm = rgm
        
        if app is not None:
            self.app = app
            
        self.view.setUrl(QUrl("params.html"))

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
    
    dsinfo = attrdict.AttrDict(
        filename="",
        subsampling = (2,2,1),
        spacing = (1,1,1),
        smethod = "Hvol",
        psmethod="HvolN",
        shape = (0,0,0),
        T=0.001,
        N=1000,
        )
    
    
    rgm = None
    dataset = None
    
    
    @pyqtSlot(result=str)
    def ds_json(self):
        return json.dumps(self.dsinfo)
    
    
    @pyqtSlot()
    def save(self):        
        pickle.dump([self.dsinfo,self.dataset,self.rgm.rg],
                    open(self.dsinfo.filename+".pickle","wb"),
                    protocol=2)       

        
    
    @pyqtSlot(str)
    def reloadModel(self,ds_json):
        
        jsinfo  = attrdict.AttrDict(**json.loads(str(ds_json)))
        
        if not jsinfo.recompute:            
            self.dsinfo,self.dataset,rg = pickle.load(open(self.dsinfo.filename+".pickle"))
            self.rgm = ReebgraphModel(self.dataset,rg,self)
        else:            
            dsinfo  = self.dsinfo        
            dsinfo.update((k, jsinfo[k]) for k in set(dsinfo).intersection(jsinfo))
        
            if isinstance(dsinfo.subsampling,unicode) :
                dsinfo.subsampling = map(int,dsinfo.subsampling.split(","))            

            if isinstance(dsinfo.spacing,unicode):
                dsinfo.spacing = map(float,dsinfo.spacing.split(","))
            
            
            self.dataset = None        
            if dsinfo.filename.endswith(".vti"):            
                vd = read_vti(self.dsinfo.filename)
                self.dataset = vd.arr[::dsinfo.subsampling[2],::dsinfo.subsampling[1],::dsinfo.subsampling[0]].copy()
                dsinfo.spacing = vd.spacing
                dsinfo.shape = vd.shape[::-1]            
                #np.savez(self.filename.replace(".vti",".npz"),self.dataset)
            
            #if self.filename.endswith(".npz"):
                #self.dataset = np.load(self.filename)["arr_0"]

        
            self.rgm = None
            if self.dataset is not None:
                
                rg = pyrg.ContourTree()
                rg.computeGrid3D(self.dataset,smethod=str(dsinfo.psmethod),N=int(dsinfo.N),T=float(dsinfo.T))
                rg.computeFeatureHierarchy(smethod=str(dsinfo.smethod))
                self.rgm = ReebgraphModel(self.dataset,rg,self)
            
            
        if self.dataset is not None:
            
            ## Create Volume Render pipeline
            self.vrPipe = VolumeRenderPipeine(self.dataset,self.vtkWidget_vr.GetRenderWindow())
            
            self.vrPipe.setSpacing(tuple(s*m for s,m in zip(self.dsinfo.spacing,self.dsinfo.subsampling)))
            
            self.rgm.dsChanged.connect(self.vrPipe.reloadData)
        
        
        self.webview.refresh(self.rgm,app=self)

        
    
    
 
    def __init__(self, parent = None, filename=""):
        QtGui.QMainWindow.__init__(self, parent)   
        
        self.dsinfo.filename = filename

                        
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
        
        #self.reloadModel("")
        self.webview.refresh(self.rgm,app=self)
        

         

def main():

    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow(filename="" if sys.argv <=1 else sys.argv[1])
 
    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()
