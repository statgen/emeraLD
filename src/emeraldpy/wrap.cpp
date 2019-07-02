#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;

#include <vector>
#include "emeraldpy/core.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<LDResult>);

PYBIND11_MODULE(emeraldpy, m) {
  m.def("ld_refsnp", &ld_refsnp, "Calculate LD for a reference variant from a m3vcf or VCF file.");
  py::class_<LDResult>(m, "LDResult")
    .def_readwrite("chrom", &LDResult::chrom)
    .def_readwrite("pos1", &LDResult::pos1)
    .def_readwrite("rsid1", &LDResult::rsid1)
    .def_readwrite("ref1", &LDResult::ref1)
    .def_readwrite("alt1", &LDResult::alt1)
    .def_readwrite("epacts1", &LDResult::epacts1)
    .def_readwrite("pos2", &LDResult::pos2)
    .def_readwrite("rsid2", &LDResult::rsid2)
    .def_readwrite("ref2", &LDResult::ref2)
    .def_readwrite("alt2", &LDResult::alt2)
    .def_readwrite("epacts2", &LDResult::epacts2)
    .def_readwrite("r", &LDResult::r)
    .def_readwrite("rsq", &LDResult::rsq)
    .def_readwrite("d", &LDResult::d)
    .def_readwrite("dprime", &LDResult::dprime);
  py::bind_vector<std::vector<LDResult>>(m, "LDResultVector");
}
