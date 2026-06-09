#define CONFORMER_TOOLKIT_NO_MAIN
#include "conformer_deduplicate.cpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(conformer_toolkit_cpp, m) {
    m.doc() = "pybind11 interface for Conformer_Batch C++ conformer deduplication";

    py::class_<Vec3>(m, "Vec3")
        .def(py::init<>())
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z);

    py::class_<Molecule>(m, "Molecule")
        .def(py::init<>())
        .def_readwrite("symbols", &Molecule::symbols)
        .def_readwrite("coords", &Molecule::coords);

    py::class_<Conformer_Record>(m, "Conformer_Record")
        .def_readonly("input_index", &Conformer_Record::input_index)
        .def_readonly("source", &Conformer_Record::source)
        .def_readonly("comment", &Conformer_Record::comment)
        .def_readonly("mol", &Conformer_Record::mol);

    py::class_<UniqueEntry>(m, "UniqueEntry")
        .def_readonly("input_index", &UniqueEntry::input_index)
        .def_readonly("path", &UniqueEntry::path)
        .def_readonly("mol", &UniqueEntry::mol);

    py::class_<DuplicateEntry>(m, "DuplicateEntry")
        .def_readonly("input_index", &DuplicateEntry::input_index)
        .def_readonly("path", &DuplicateEntry::path)
        .def_readonly("representative_input_index", &DuplicateEntry::representative_input_index)
        .def_readonly("representative_path", &DuplicateEntry::representative_path)
        .def_readonly("rmsd", &DuplicateEntry::rmsd);

    py::class_<Conformer_Batch::Remove_Duplicates_Result>(m, "Remove_Duplicates_Result")
        .def_readonly("unique", &Conformer_Batch::Remove_Duplicates_Result::unique)
        .def_readonly("duplicates", &Conformer_Batch::Remove_Duplicates_Result::duplicates);

    py::class_<Conformer_Batch>(m, "Conformer_Batch")
        .def_static("from_xyz_files", &Conformer_Batch::from_xyz_files, py::arg("paths"))
        .def_static(
            "from_xyz_directory",
            [](const std::string& directory) {
                return Conformer_Batch::from_xyz_directory(directory);
            },
            py::arg("directory")
        )
        .def_static(
            "from_multi_xyz",
            [](const std::string& path) {
                return Conformer_Batch::from_multi_xyz(path);
            },
            py::arg("path")
        )
        .def("size", &Conformer_Batch::size)
        .def("records", &Conformer_Batch::records, py::return_value_policy::reference_internal)
        .def(
            "index_cleanup",
            &Conformer_Batch::index_cleanup,
            py::arg("max_mappings") = 1000000,
            py::arg("bond_scale") = 1.1,
            py::arg("ambiguity_gap") = 1e-6
        )
        .def(
            "remove_duplicates",
            &Conformer_Batch::remove_duplicates,
            py::arg("tolerance") = 1e-3,
            py::arg("run_index_cleanup") = false,
            py::arg("max_mappings") = 1000000,
            py::arg("bond_scale") = 1.1,
            py::arg("ambiguity_gap") = 1e-6
        )
        .def(
            "write_records",
            [](const Conformer_Batch& batch, const std::string& outdir, const std::string& prefix) {
                batch.write_records(outdir, prefix);
            },
            py::arg("outdir"),
            py::arg("prefix") = "cleaned"
        );
}
