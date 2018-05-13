#include <iostream>
#include <sstream>


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "utl.h"
#include "Grid3D.hpp"
#include "MergeTree.hpp"
#include "ContourTreeData.hpp"
#include "ContourTreeData.hpp"
#include "SimplifyCT.hpp"
#include "Persistence.hpp"

using namespace contourtree;


namespace py=pybind11;

py::tuple simplifyContourTree(ContourTree &ct) {

    std::vector<int64_t> nodeids;
    std::vector<scalar_t> nodefns;
    std::vector<char> nodeTypes;
    std::vector<int64_t> arcs;
    std::vector<uint32_t> arcMap;
	arcMap.resize(ct.nv);

	ct.output(nodeids,nodefns,nodeTypes,arcs,arcMap.data());

    std::cout << SVAR(nodeids.size()) << std::endl;
    std::cout << SVAR(nodefns.size()) << std::endl;
    std::cout << SVAR(nodeTypes.size()) << std::endl;
    std::cout << SVAR(arcs.size()) << std::endl;
    std::cout << SVAR(arcMap.size()) << std::endl;


    ContourTreeData ctdata;
    ctdata.loadData(nodeids,nodefns,nodeTypes,arcs);

    SimplifyCT sim;
    sim.setInput(&ctdata);

    Persistence simfn(ctdata);
    sim.simplify(&simfn);

    std::vector<uint32_t> order;
    std::vector<float>   wts;

    sim.outputOrder(order,wts);

    return py::make_tuple(py::array(order.size(),order.data()),
                          py::array(wts.size(),wts.data()));



}


struct py_icoord {
    int x;
    int y;
    int z;
};

struct py_node {
    int64_t    id;
    py_icoord  crd;
    scalar_t   fn;
    char       type;
};


py::tuple computeCT_Grid3D(py::array_t<scalar_t> grid){

    size_t X = grid.shape(2);
    size_t Y = grid.shape(1);
    size_t Z = grid.shape(0);


	Grid3D g(grid.mutable_data(),X,Y,Z);

    MergeTree mt;
    mt.computeTree(&g,TypeContourTree);

    std::vector<int64_t> nodeids;
    std::vector<scalar_t> nodefns;
    std::vector<char> nodeTypes;
    std::vector<int64_t> arcs;
	py::array_t<uint32_t> arcMap;
	arcMap.resize({Z,Y,X});

	mt.ctree.output(nodeids,nodefns,nodeTypes,arcs,arcMap.mutable_data());

    std::vector<py_node> nodes(nodeids.size());

    std::map<int64_t,int64_t> idToNodeNum;

    for(int i = 0; i < nodes.size(); ++i ) {
        int id = nodeids[i];
        nodes[i].id = id;
        nodes[i].crd.x = id%X;
        nodes[i].crd.y = (id/X)%Y;
        nodes[i].crd.z = (id/(X*Y))%Z;
        nodes[i].fn = nodefns[i];
        nodes[i].type = nodeTypes[i];
        idToNodeNum[id] = i;
    }

    for(auto &a: arcs) {
        a = idToNodeNum[a];
    }

    std::cout << SVAR(arcMap.size()) << std::endl;

    auto nodes_   = py::array(nodes.size(),nodes.data());
    auto arcs_    = py::array(arcs.size(),arcs.data());
    arcs_.resize({arcs.size()/2,size_t(2)},true);


	return py::make_tuple(nodes_,arcs_,arcMap);
}


PYBIND11_MODULE(pyrg, m) {

    m.doc() = "Py ReebGraph module"; // optional module docstring


    PYBIND11_NUMPY_DTYPE(py_icoord,x,y,z);
    PYBIND11_NUMPY_DTYPE(py_node,id,crd,fn,type);



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
            .def("computeContourTree",[](MergeTree& mt,Grid3D& grid)->void
            {std::cout << SVAR(&mt) << std::endl;mt.computeTree(&grid,TypeContourTree);})

            .def("computeJoinTree",[](MergeTree& ct,Grid3D& grid)->void{ct.computeTree(&grid,TypeJoinTree);})
            .def("computeSplitTree",[](MergeTree& mt,Grid3D& grid)->void{
                std::cout << SVAR(&mt) << std::endl;
                mt.computeTree(&grid,TypeSplitTree);})
            .def("genPersistenceHierarchy",[](MergeTree &mt){
                std::cout << SVAR(&mt) << std::endl;
                return simplifyContourTree(mt.ctree);})
    ;

    m.def("computeCT_Grid3D",&computeCT_Grid3D,"Compute the contour tree on a structured3D grid");


}
