#include <iostream>
#include <sstream>

#include <pybind11/pybind11.h>

#include "../utl.h"
#include "../DisjointSets.hpp"

#include "test.hpp"

using namespace contourtree;
namespace py=pybind11;


PYBIND11_MODULE(pyrgtest, m) {
    m.doc() = "Py ReebGraph module"; // optional module docstring

    m.def("testDisjointSets", &testDisjointSets, "testDisjointSets");
}
