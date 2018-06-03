import os,sys

from PyQt4.QtCore import *
from PyQt4 import  QtGui
from PyQt4.QtWebKit import *

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import pyrg
import vtk

from ex2 import create_3gauss,VolumeRenderPipeine,read_vti


class ReebgraphModel(QObject):
    """
    Class that handles creating and maintaining the reebgraph
    Also provides jsonified data for webview
    """
    
    dsChanged = pyqtSignal()    

    def __init__(self, dataset, parent = None):
        QObject.__init__(self, parent)
        self.ods = np.copy(dataset)
        self.ds  = dataset
        
        self.rg = pyrg.ContourTree()
        self.rg.computeGrid3D(self.ds)
        self.rg.computeFeatureHierarchy()
                
        print len(self.rg.nodes)
        print len(self.rg.arcs)
        
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
        import json
        
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
        import json
        
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
        
        
    def refresh(self,rgm=None):
        
        if rgm != None:
            self.rgm = rgm        
                        
        self.view.setUrl(QUrl("force2.html"))

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
    
 
    def __init__(self, parent = None, dataset=None):
        QtGui.QMainWindow.__init__(self, parent)
        
        # Create Dataset
        self.dataset = create_3gauss() if dataset is None else dataset   
                        
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
        
        self.rgm = ReebgraphModel(self.dataset,self)
        
        self.webview.refresh(self.rgm)
        
        ## Create Volume Render pipeline
        self.volumeRenderPipeline = VolumeRenderPipeine(self.dataset,self.vtkWidget_vr.GetRenderWindow())
        
        self.rgm.dsChanged.connect(self.volumeRenderPipeline.reloadData)

         

def main():
    dataset = None
    
    if len(sys.argv) >1:
        fn = sys.argv[1]
        
        if fn.endswith(".vti"):
            dataset = read_vti(fn)
            np.savez(fn.replace(".vti",".npz"),dataset)
            
        if fn.endswith(".npz"):
            dataset = np.load(fn)["arr_0"]

    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow(dataset=dataset)
 
    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()
