import os,sys

from PyQt4.QtCore import *
from PyQt4 import  QtGui
from PyQt4.QtWebKit import *

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import pyrg
import vtk

from ex2 import create_3gauss,VolumeRenderPipeine


class ReebgraphModel(QObject):
    """
    Class that handles creating and maintaining the reebgraph
    Also provides jsonified data for webview
    """

    def __init__(self, dataset, parent = None):
        QObject.__init__(self, parent)
        self.ds  = dataset        
        self.nodes,self.arcs,self.arcmap  = pyrg.computeCT_Grid3D(self.ds)
        self.sorder,self.swts,self.shier  = pyrg.simplifyCT_Pers(self.nodes,self.arcs)           
    
    @pyqtSlot(result=str)
    def json(self):
        import json
        
        nodes,arcs = self.nodes,self.arcs

        rng   = [float(nodes["fn"].min()),float(nodes["fn"].max())]                
        nodes = [ {"id":i, "name":str(n["id"]), "fn":float(n["fn"]),"group":int(n["type"]) } for i,n in enumerate(nodes)]
        links = [ {"source":arc[0], "target":arc[1]} for arc in arcs ]
        
        for (a,b),w in zip(self.sorder,self.swts):
            nodes[a]["weight"] = float(w);
            nodes[b]["weight"] = float(w);

        
        return json.dumps({"nodes":nodes,"links":links,"range":rng})

    @pyqtSlot(result=str)
    def hierTree_json(self):
        import json
        
        nodes,arcs = self.nodes,self.arcs
        
        def get_type(a,b):
            return ("min" if nodes[a]["type"] == 1 else "sad") +"-" + ("max" if nodes[b]["type"] == 2 else "sad")
        
        arcTree = dict([(str((a,b)),
                         {"name":str((a,b)),
                          #"size":float(nodes[b]["fn"] - nodes[a]["fn"]),
                          #"size":10,
                          "size":1 + 100*float(nodes[b]["fn"] - nodes[a]["fn"]),
                          "type":get_type(a,b),
                          "weight":1.0
                          }) for (a,b) in arcs])
                                 
        for (c,m,d,u),w in zip(self.shier,self.swts):
                        
            clu = str((min(c,m),max(c,m)))
            m_u = str((m,u))
            d_m = str((d,m))
            d_u = str((d,u))
            
            arcTree[clu]["weight"] = arcTree[m_u]["weight"] = arcTree[d_m]["weight"] = float(w)
            
            assert not arcTree.has_key((d_u))
            arcTree[d_u] = {
                "name":d_u,
                "children": [arcTree[clu],arcTree[d_m],arcTree[m_u]],
                #"size":float(nodes[u]["fn"] - nodes[d]["fn"]),
                "type":get_type(d,u),
                "weight":1.0,                
                }
            del arcTree[m_u]
            del arcTree[d_m]
            del arcTree[clu]
            
        
        for k,v in arcTree.iteritems():
            arcTree = {"name":str(k),
                       "children":v["children"],
                       "type":v["type"]
                       }
            break
        
        return json.dumps(arcTree)

    
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
                        
        self.view.setUrl(QUrl("packLayout.html"))

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
        self.dataset = create_3gauss();   
                        
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
        
        self.webview.refresh(ReebgraphModel(self.dataset,self))
        
        ## Create Volume Render pipeline
        self.volumeRenderPipeline = VolumeRenderPipeine(self.dataset,self.vtkWidget_vr.GetRenderWindow()) 
         
        


def main():
    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow()
 
    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()
