# ReebGraph
ReebGraph implementation based On HarishD's Code/Algorithms

The main aim is to provide a good python interface for efficient computation of ReebGraphs of 3d structured grids


## Installation

The installation is tested on Ubuntu 14.04 and 16.04. 

You will need python2.7 with python-dev header files installed. 
Cmake > 3.9 (Ubuntu 16.04 comes with 3.5. you can download 3.9 from [here](https://cmake.org/download/) ). 

Install the code as follows. 

```bash
$ git clone https://github.com/nithins/ReebGraph.git
$ cd ReebGraph
$ git submodule update --init --recursive
$ mkdir build
$ cd build
$ ~/software/cmake-3.10.2-Linux-x86_64/bin/cmake ../
$ make -j8
```
Invariabley, cmake will not find the python dev paths. Run (the respective 3.9+) cmake-gui and adjust the paths. 

## Running the Code

After the build, a file name pyrg.so will be placed in the example folder to run the examples. They must be run from the example folder itself. 
```bash
$ cd ../examples
```

### Example 1

This computes the Contour tree of a simple synthetic function (3 3D gaussians) in a unit volume

```bash
$ python ex1.py
```
### Example 2

This renders the Contour Tree along with a volume rendering vtk and qt4. (Contour tree rendering is incomplete right now)

```bash
$ python ex2.py
```

This can also compute the CT of a given .vti file. If this is not given, the default 3gauss dataset of ex1 is used. 

```bash
$ python ex2.py file.vti
```


Note: this requires a working installation of pyqt with vtk built with python and qt support. When last tested, Ubuntu's 16.04's default installation of vtk python packages did not work properly. 

VTK 8.0 worked by compiling from source.  The qt version was set to 5 (even so it builds the qt4 viewer widget) and python version to 2.7. Ubuntu packages for qt-dev as well as pyqt4 were installed. 


### Example 3

This Example plots the CT using various layouts from [d3.js](https://d3js.org/) . 
Currently, the packLayout/circleLayout is rendered to show the feature tree. 
It can be used to pick components and see which parts of the Volume they corresspond to. 
Simplification thresholds may be used to interactively show/hide nodes of the feature tree. 

As with the previous example, it may be used to compute for a dataset given int .vti format. (or the default 3gauss if none is given). 

```bash
$ python ex3.py file.vti
```



