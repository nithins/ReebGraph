#include <iostream>
#include <sstream>

#include <pybind11/pybind11.h>

#include "test.hpp"

using namespace contourtree;
namespace py=pybind11;


PYBIND11_MODULE(pyrgtest, m) {
    m.doc() = "Py ReebGraph module"; // optional module docstring

    m.def("testDisjointSets", &testDisjointSets, "testDisjointSets");
    m.def("testToyDataset", []{
        generateData("toy.raw");
        toyProcessing("toy.raw");
        toyFeatures("toy.raw");
    }, "testToyDataset");
}
