#include <iostream>
#include <sstream>


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "utl.h"
#include "Grid3D.hpp"
#include "MergeTree.hpp"
#include "ContourTreeData.hpp"
#include "ContourTreeData.hpp"
#include "SimplifyCT.hpp"
#include "SimplifyCT2.hpp"
#include "Persistence.hpp"

using namespace contourtree;
using namespace pybind11::literals;

namespace py=pybind11;

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



typedef contourtree::arc_t arc_t;

typedef contourtree::arc_t py_arc;

struct pyContourTree {
	py::array_t<py_node>  nodes;
	py::array_t<int64_t>  arcs;
	py::array_t<uint32_t> arcMap;


	py::array_t<float>    fwts;
	py::array_t<uint32_t> fhier;
	py::array_t<int64_t>  farcs;


    static py::tuple pickle(const pyContourTree & pct){
        return py::make_tuple(
                    pct.nodes,
                    pct.arcs,
                    pct.arcMap,
                    pct.fwts,
                    pct.fhier,
                    pct.farcs);
    }

    static pyContourTree unpickle(py::tuple t){
        ENSURES(t.size() == 6);
        pyContourTree pct;
        pct.nodes  = t[0].cast<decltype(nodes)>();
        pct.arcs   = t[1].cast<decltype(arcs)>();
        pct.arcMap = t[2].cast<decltype(arcMap)>();
        pct.fwts   = t[3].cast<decltype(fwts)>();
        pct.fhier  = t[4].cast<decltype(fhier)>();
        pct.farcs  = t[5].cast<decltype(farcs)>();
        return pct;
    }



    bool computeGrid3D(py::array_t<scalar_t> &grid,std::string smethod,float T, int N)
    {

        ENSURES(smethod == "PersT" || smethod == "PersN" || smethod == "HvolN" || smethod == "" ) << SVAR(smethod);


		size_t X = grid.shape(2);
		size_t Y = grid.shape(1);
		size_t Z = grid.shape(0);


		Grid3D g(grid.mutable_data(),X,Y,Z);

		MergeTree mt;
		mt.computeTree(&g,TypeContourTree);

		std::vector<int64_t> nodeids;
		std::vector<scalar_t> nodefns;
		std::vector<char> nodeTypes;


		std::vector<py_node> nodes;
		std::vector<arc_t> arcs;


		arcMap.resize({Z,Y,X});
		auto arcMapPtr = arcMap.mutable_data();

		{
			std::vector<int64_t> carcs;

			mt.ctree.output(nodeids,nodefns,nodeTypes,carcs,arcMapPtr);

			std::map<int64_t,int64_t> idToNodeNum;
			for(int i = 0; i < nodeids.size(); ++i )
				idToNodeNum[nodeids[i]] = i;

			for(int i = 0 ; i < carcs.size(); i+=2)
				arcs.push_back({idToNodeNum[carcs[i]],idToNodeNum[carcs[i+1]]});
		}


        auto arcRemap = contourtree::splitMonkeysAndNazis(nodeids,nodefns,nodeTypes,arcs);

        if(!smethod.empty()) {

            std::vector<float> arcVols(arcs.size(),0);

            if(smethod == "HvolN") {
                for(int i = 0 ; i < arcMap.size(); ++i){
                    auto ai = *(arcMap.data() + i);
                    ENSURES( 0 <= ai && ai < arcRemap.size());
                    arcVols[arcRemap[ai]]++;
                }
            }

            auto arcRemap2 = contourtree::preSimplify(nodeids,nodefns,nodeTypes,arcs,arcVols,smethod,T,N);

            for(int i = 0 ; i < arcRemap.size(); ++i)
                arcRemap[i] = arcRemap2[arcRemap[i]];

		}


        for(int i = 0 ; i < X*Y*Z; ++i)
            arcMapPtr[i] = arcRemap[arcMapPtr[i]];



		nodes.resize(nodeids.size());

		for(int i = 0; i < nodes.size(); ++i ) {
			int id = nodeids[i];
			nodes[i].id = id;
			nodes[i].crd.x = id%X;
			nodes[i].crd.y = (id/X)%Y;
			nodes[i].crd.z = (id/(X*Y))%Z;
			nodes[i].fn = nodefns[i];
			nodes[i].type = nodeTypes[i];
		}

		this->nodes = py::array(nodes.size(),nodes.data());
		this->arcs  = py::array({arcs.size(),size_t(2)},&arcs[0].first);


		return true;
	}

    bool computeFeatureHierarchy(std::string smethod) {

		std::vector<int64_t> nodeids;
		std::vector<scalar_t> nodefns;
		std::vector<char> nodeTypes;
		std::vector<arc_t> arcs;

		for(int i = 0; i < nodes.size(); ++i ) {
			nodeids.push_back(nodes.at(i).id);
			nodefns.push_back(nodes.at(i).fn);
			nodeTypes.push_back(nodes.at(i).type);
		}

		for(int i = 0 ; i < this->arcs.shape(0); ++i) {
			arcs.push_back({this->arcs.at(i,0),this->arcs.at(i,1)});
		}

		std::vector<float>    wts;
		std::vector<arc_t>    carcs;
		std::vector<uint32_t> featureHierarchy;
		std::vector<arc_t>    sarcs;


        if(smethod == "Pers") {


			contourtree::simplifyPers(nodefns,nodeTypes,arcs,carcs,wts,featureHierarchy,sarcs);

		}
        else if(smethod == "Hvol")
        {
			std::vector<float> arcVols(arcs.size(),0);

			auto numArcMap = arcMap.size();

			for(int i = 0 ; i < numArcMap; ++i){
				auto ai = *(arcMap.data() + i);
				ENSURES( 0 <= ai && ai < arcs.size());
				arcVols[ai]++;
			}

			contourtree::simplifyHyperVolume(nodefns,nodeTypes,arcs,arcVols,
											 carcs,wts,featureHierarchy,sarcs);

		}
        else {
            ENSURES(false) << "Unknown simp method" << SVAR(smethod) ;
        }

		ENSURES(sarcs.size() == 1) <<SVAR(sarcs.size());
		ENSURES(wts.size() == carcs.size() && wts.size() == (featureHierarchy.size()/5))
				<<SVAR(wts.size()) << SVAR(carcs.size()) << SVAR(featureHierarchy.size());

		wts.push_back(1.0);
		carcs.push_back({sarcs.front().first,sarcs.back().second});


		this->farcs = py::array({wts.size(),(size_t)2},&carcs[0].first);
		this->fwts  = py::array(wts.size(),wts.data());
		this->fhier = py::array({wts.size()-1,(size_t)5},featureHierarchy.data());
		return true;
	}
};



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

    py::class_<pyContourTree>(m, "ContourTree","Computation of Contour trees and their features")
			.def(py::init<>())
            .def("computeGrid3D",&pyContourTree::computeGrid3D,"grid"_a,"smethod"_a="PersT","T"_a=0.001,"N"_a=1000)
            .def_readonly("nodes",&pyContourTree::nodes,"nodes of the contour tree")
            .def_readonly("arcs",&pyContourTree::arcs,"arcs of the contour tree")
            .def_readonly("arcmap",&pyContourTree::arcMap,"mapping from vertices of the domain into arcs ")

            .def("computeFeatureHierarchy",&pyContourTree::computeFeatureHierarchy,"smethod"_a="Hvol")
            .def_readonly("farcs",&pyContourTree::farcs,"Feature arcs a.k.a brach decomposition")
            .def_readonly("fhier",&pyContourTree::fhier,"Feature hierarchy ")
            .def_readonly("fwts",&pyContourTree::fwts,"Feature Weights ")

            .def(py::pickle(&pyContourTree::pickle,&pyContourTree::unpickle))
            ;


}
