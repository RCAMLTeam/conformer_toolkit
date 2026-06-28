#define CONFORMER_TOOLKIT_NO_MAIN
#include "conformer_deduplicate.cpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_native, m) {
    m.doc() = "pybind11 interface for C++ conformer deduplication and molecule ring analysis";

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
        .def_readonly("properties", &Conformer_Record::properties)
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

    py::class_<Ring_Atom_Adjacency>(
        m,
        "Ring_Atom_Adjacency",
        "Adjacency for one ring atom. adjacent_atoms contains ring neighbors and directly bonded non-ring atoms."
    )
        .def_readonly("atom", &Ring_Atom_Adjacency::atom)
        .def_readonly("adjacent_atoms", &Ring_Atom_Adjacency::adjacent_atoms);

    py::class_<Ring_Record>(
        m,
        "Ring_Record",
        "One detected ring, with spatially ordered atom indices and per-atom molecule adjacency."
    )
        .def_readonly("atoms", &Ring_Record::atoms)
        .def_readonly("adjacency_list", &Ring_Record::adjacency_list);

    py::class_<Conformer_Group::Remove_Duplicates_Result>(m, "Remove_Duplicates_Result")
        .def_readonly("unique", &Conformer_Group::Remove_Duplicates_Result::unique)
        .def_readonly("duplicates", &Conformer_Group::Remove_Duplicates_Result::duplicates);

    py::class_<Conformer_Group>(
        m,
        "Conformer_Group",
        "Container for conformers of one molecule, with deduplication and RDKit-derived ring information."
    )
        .def_static(
            "from_xyz_files",
            &Conformer_Group::from_xyz_files,
            py::arg("paths"),
            py::arg("comment_template") = "",
            "Load one conformer from each XYZ file path."
        )
        .def_static(
            "from_xyz_directory",
            [](const std::string& directory, const std::string& comment_template) {
                return Conformer_Group::from_xyz_directory(directory, comment_template);
            },
            py::arg("directory"),
            py::arg("comment_template") = "",
            "Load all .xyz files in a directory, sorted by path."
        )
        .def_static(
            "from_multi_xyz",
            [](const std::string& path, const std::string& comment_template) {
                return Conformer_Group::from_multi_xyz(path, comment_template);
            },
            py::arg("path"),
            py::arg("comment_template") = "",
            "Load conformers from a concatenated multi-XYZ file."
        )
        .def("size", &Conformer_Group::size)
        .def("__len__", &Conformer_Group::size)
        .def("sort_by_energy", &Conformer_Group::sort_by_energy,
            py::arg("energy_property") = "energy",
            "Sort conformers in place from lowest to highest energy (stable for ties).")
        .def("filter_by_maximum_energy", &Conformer_Group::filter_by_maximum_energy,
            py::arg("maximum_energy"), py::arg("energy_property") = "energy",
            "Keep conformers whose energy is less than or equal to maximum_energy.")
        .def("retain_lowest_energy_percent", &Conformer_Group::retain_lowest_energy_percent,
            py::arg("percent"), py::arg("energy_property") = "energy",
            "Sort by energy and retain ceil(size * percent / 100) conformers.")
        .def("filter_by_boltzmann_population_ratio", &Conformer_Group::filter_by_boltzmann_population_ratio,
            py::arg("minimum_ratio"), py::arg("temperature_kelvin") = 298.15,
            py::arg("energy_property") = "energy", py::arg("energy_to_joules_per_mole") = 1000.0,
            "Keep conformers with Boltzmann population at least minimum_ratio times that of the energy minimum.")
        .def(
            "records",
            &Conformer_Group::records,
            py::return_value_policy::reference_internal,
            "Return stored conformer records."
        )
        .def(
            "index_cleanup",
            &Conformer_Group::index_cleanup,
            py::arg("max_mappings") = 1000000,
            py::arg("bond_scale") = 1.3,
            py::arg("ambiguity_gap") = 1e-6,
            py::arg("charge") = 0,
            "Normalize all conformers to the atom indexing of the first conformer."
        )
        .def(
            "remove_duplicates",
            &Conformer_Group::remove_duplicates,
            py::arg("tolerance") = 1e-3,
            py::arg("run_index_cleanup") = false,
            py::arg("max_mappings") = 1000000,
            py::arg("bond_scale") = 1.3,
            py::arg("ambiguity_gap") = 1e-6,
            py::arg("charge") = 0,
            "Remove duplicate conformers by chemistry-aware RMSD comparison."
        )
        .def(
            "write_records",
            [](const Conformer_Group& group, const std::string& outdir, const std::string& prefix) {
                group.write_records(outdir, prefix);
            },
            py::arg("outdir"),
            py::arg("prefix") = "cleaned",
            "Write the currently stored conformer records as XYZ files."
        )
        .def(
            "rings",
            &Conformer_Group::rings,
            py::return_value_policy::reference_internal,
            "Return ring records populated by detect_rings()."
        )
        .def(
            "adjacency_table",
            &Conformer_Group::adjacency_table,
            py::return_value_policy::reference_internal,
            "Return the molecule adjacency table inferred from conformer 0 by detect_rings()."
        )
        .def(
            "detect_rings",
            &Conformer_Group::detect_rings,
            py::arg("bond_scale") = 1.3,
            py::arg("charge") = 0,
            "Infer chemistry with RDKit from the first conformer and store adjacency/ring information."
        );

    m.attr("Conformer_Batch") = m.attr("Conformer_Group");
}
