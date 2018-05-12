#include <iostream>
#include <sstream>


#include <pybind11/pybind11.h>

//using namespace contourtree;


namespace py=pybind11;


PYBIND11_MODULE(pyrg, m) {
    m.doc() = "Py ReebGraph module"; // optional module docstring

    // m.def("testDisjointSets", &testDisjointSets, "testDisjointSets");
}
