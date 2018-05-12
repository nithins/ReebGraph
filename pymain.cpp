#include <iostream>
#include <sstream>


#include <pybind11/pybind11.h>

#include "utl.h"
#include "Grid3D.hpp"
#include "MergeTree.hpp"

using namespace contourtree;


namespace py=pybind11;


PYBIND11_MODULE(pyrg, m) {

    m.doc() = "Py ReebGraph module"; // optional module docstring



    py::class_<Grid3D>(m, "Grid3D", py::buffer_protocol())
            .def(py::init([](py::tuple shape){ENSURES(shape.size() == 3);
                return new Grid3D(shape[2].cast<int>(),shape[1].cast<int>(),shape[0].cast<int>());}))
            .def("atXYZ",[](Grid3D &m, int x,int y, int z){return m.getFunctionValue(m.index(x,y,z));})
            .def_buffer([](Grid3D &m) -> py::buffer_info {
                return py::buffer_info(
                    m.data(),                                  /* Pointer to buffer */
                    sizeof(scalar_t),                          /* Size of one scalar */
                    py::format_descriptor<scalar_t>::format(), /* Python struct-style format descriptor */
                    3,                                         /* Number of dimensions */
                    { m.dimZ(), m.dimY(),m.dimX()},            /* Buffer dimensions */
                    { sizeof(scalar_t) * m.dimX()*m.dimY(),    /* Strides (in bytes) for each index */
                      sizeof(scalar_t) * m.dimX(),
                      sizeof(scalar_t) });
            })
    ;


    py::class_<MergeTree>(m, "MergeTree", py::buffer_protocol())
            .def(py::init<>())
            .def("computeContourTree",[](MergeTree& ct,Grid3D& grid)->void{ct.computeTree(&grid,TypeContourTree);})
            .def("computeJoinTree",[](MergeTree& ct,Grid3D& grid)->void{ct.computeTree(&grid,TypeJoinTree);})
            .def("computeSplitTree",[](MergeTree& ct,Grid3D& grid)->void{ct.computeTree(&grid,TypeSplitTree);})
    ;


}
