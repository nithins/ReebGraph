#include <iostream>
#include <sstream>

#include <pybind11/pybind11.h>

#include "../utl.h"
#include "../DisjointSets.hpp"

using namespace contourtree;
namespace py=pybind11;


void testDisjointSets() {
    int numElements = 128;
    int numInSameSet = 16;
    DisjointSets<int64_t> ds(numElements);
    int set1, set2;

    for (int k = 1; k < numInSameSet; k *= 2) {
        for (int j = 0; j + k < numElements; j += 2 * k) {
            set1 = ds.find(j);
            set2 = ds.find(j + k);
            ds.merge(set1, set2);
        }
    }

    for (int i = 0; i < numElements; i++) {
        std::cout << ds.find(i) << "*";
        if (i % numInSameSet == numInSameSet - 1)
            std::cout << "\n";

        ENSURES(ds.find(i) == numInSameSet*int(i/numInSameSet))
                << " ds.find(i)=" << ds.find(i)
                << " int(i/numInSameSet)" << int(i/numInSameSet);
    }
    std::cout << "\n";
}



PYBIND11_MODULE(pyrgtest, m) {
    m.doc() = "Py ReebGraph module"; // optional module docstring

    m.def("testDisjointSets", &testDisjointSets, "testDisjointSets");
}
