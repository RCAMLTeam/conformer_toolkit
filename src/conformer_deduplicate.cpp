#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Vec3 {
    double x{};
    double y{};
    double z{};
};

struct Molecule {
    std::vector<std::string> symbols;
    std::vector<Vec3> coords;
};

struct Conformer_Record {
    std::size_t input_index{};
    std::string source;
    std::string comment;
    Molecule mol;
};

struct MolecularGraph {
    std::vector<std::vector<bool>> bonded;
    std::vector<std::string> signatures;
};

struct UniqueEntry {
    std::size_t input_index{};
    std::string path;
    Molecule mol;
};

struct DuplicateEntry {
    std::size_t input_index{};
    std::string path;
    std::size_t representative_input_index{};
    std::string representative_path;
    double rmsd{};
};

struct RmsdMatch {
    double best{std::numeric_limits<double>::infinity()};
    double second{std::numeric_limits<double>::infinity()};
    std::size_t mappings_checked{};
    std::vector<std::size_t> best_mapping;
};

static bool read_xyz_record(std::istream& in, Molecule& mol, std::string& comment) {
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty()) {
            break;
        }
    }

    if (!in && line.empty()) {
        return false;
    }

    std::istringstream count_stream(line);
    std::size_t atom_count = 0;
    if (!(count_stream >> atom_count)) {
        throw std::runtime_error("XYZ record first line is not an atom count");
    }

    if (!std::getline(in, comment)) {
        throw std::runtime_error("XYZ record is missing comment line");
    }

    mol = Molecule{};
    mol.symbols.reserve(atom_count);
    mol.coords.reserve(atom_count);

    for (std::size_t i = 0; i < atom_count; ++i) {
        if (!std::getline(in, line)) {
            throw std::runtime_error("XYZ record ended before all atoms were read");
        }

        std::istringstream atom_stream(line);
        std::string symbol;
        Vec3 coord;
        if (!(atom_stream >> symbol >> coord.x >> coord.y >> coord.z)) {
            throw std::runtime_error("Bad XYZ atom line: " + line);
        }

        mol.symbols.push_back(symbol);
        mol.coords.push_back(coord);
    }

    return true;
}

static void write_xyz(const std::filesystem::path& path, const Molecule& mol, const std::string& comment) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not write XYZ file: " + path.string());
    }

    out << mol.coords.size() << "\n";
    out << comment << "\n";
    out << std::fixed << std::setprecision(10);
    for (std::size_t i = 0; i < mol.coords.size(); ++i) {
        out << mol.symbols[i] << " " << mol.coords[i].x << " " << mol.coords[i].y << " " << mol.coords[i].z << "\n";
    }
}

static Vec3 centroid(const std::vector<Vec3>& coords) {
    Vec3 c;
    for (const Vec3& v : coords) {
        c.x += v.x;
        c.y += v.y;
        c.z += v.z;
    }

    const double inv_n = 1.0 / static_cast<double>(coords.size());
    c.x *= inv_n;
    c.y *= inv_n;
    c.z *= inv_n;
    return c;
}

static std::vector<Vec3> centered(const std::vector<Vec3>& coords) {
    const Vec3 c = centroid(coords);
    std::vector<Vec3> out;
    out.reserve(coords.size());
    for (const Vec3& v : coords) {
        out.push_back({v.x - c.x, v.y - c.y, v.z - c.z});
    }
    return out;
}

static double norm_squared_sum(const std::vector<Vec3>& coords) {
    double total = 0.0;
    for (const Vec3& v : coords) {
        total += v.x * v.x + v.y * v.y + v.z * v.z;
    }
    return total;
}

static double covalent_radius(const std::string& symbol) {
    // Covalent radii are only used to infer a rough graph from XYZ.
    // If --allow-reorder misses expected duplicates, --bond-scale and this table
    // are the first places to inspect.
    if (symbol == "H") return 0.31;
    if (symbol == "B") return 0.84;
    if (symbol == "C") return 0.76;
    if (symbol == "N") return 0.71;
    if (symbol == "O") return 0.66;
    if (symbol == "F") return 0.57;
    if (symbol == "P") return 1.07;
    if (symbol == "S") return 1.05;
    if (symbol == "Cl") return 1.02;
    if (symbol == "Br") return 1.20;
    if (symbol == "I") return 1.39;
    if (symbol == "Si") return 1.11;
    if (symbol == "Na") return 1.66;
    if (symbol == "K") return 2.03;
    if (symbol == "Li") return 1.28;
    if (symbol == "Mg") return 1.41;
    if (symbol == "Ca") return 1.76;
    if (symbol == "Zn") return 1.22;
    if (symbol == "Fe") return 1.32;
    if (symbol == "Cu") return 1.32;
    throw std::runtime_error("Covalent radius not defined for element: " + symbol);
}

static double distance(const Vec3& a, const Vec3& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static MolecularGraph infer_graph(const Molecule& mol, double bond_scale) {
    // XYZ has no bond table, so this builds a debugging-friendly approximation:
    // distance <= bond_scale * (r_cov_a + r_cov_b) means bonded.
    const std::size_t n = mol.coords.size();
    MolecularGraph graph;
    graph.bonded.assign(n, std::vector<bool>(n, false));
    graph.signatures.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double cutoff = bond_scale * (covalent_radius(mol.symbols[i]) + covalent_radius(mol.symbols[j]));
            if (distance(mol.coords[i], mol.coords[j]) <= cutoff) {
                graph.bonded[i][j] = true;
                graph.bonded[j][i] = true;
            }
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        std::vector<std::string> neighbor_symbols;
        for (std::size_t j = 0; j < n; ++j) {
            if (graph.bonded[i][j]) {
                neighbor_symbols.push_back(mol.symbols[j]);
            }
        }
        std::sort(neighbor_symbols.begin(), neighbor_symbols.end());

        // The signature narrows atom mapping candidates before recursion.
        // It intentionally stays local and simple; RDKit mode is preferred when
        // bond order, aromaticity, charge, or stereochemistry matters.
        std::ostringstream signature;
        signature << mol.symbols[i] << "|degree=" << neighbor_symbols.size() << "|neighbors=";
        for (const std::string& neighbor : neighbor_symbols) {
            signature << neighbor << ",";
        }
        graph.signatures[i] = signature.str();
    }

    return graph;
}

static double largest_eigenvalue_symmetric_4x4(const std::array<std::array<double, 4>, 4>& matrix) {
    std::array<double, 4> q{1.0, 0.0, 0.0, 0.0};

    for (int iter = 0; iter < 100; ++iter) {
        std::array<double, 4> next{};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                next[i] += matrix[i][j] * q[j];
            }
        }

        double norm = 0.0;
        for (double value : next) {
            norm += value * value;
        }
        norm = std::sqrt(norm);
        if (norm == 0.0) {
            return 0.0;
        }

        for (double& value : next) {
            value /= norm;
        }

        double delta = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double diff = next[i] - q[i];
            delta += diff * diff;
        }
        q = next;
        if (delta < 1e-28) {
            break;
        }
    }

    double eigenvalue = 0.0;
    for (int i = 0; i < 4; ++i) {
        double row_dot = 0.0;
        for (int j = 0; j < 4; ++j) {
            row_dot += matrix[i][j] * q[j];
        }
        eigenvalue += q[i] * row_dot;
    }
    return eigenvalue;
}

static double aligned_rmsd(const Molecule& a, const Molecule& b) {
    if (a.coords.empty()) {
        return 0.0;
    }

    // Horn quaternion RMSD after removing translation. This assumes atom i in
    // both Molecule objects is already the intended mapping.
    const std::vector<Vec3> p = centered(a.coords);
    const std::vector<Vec3> q = centered(b.coords);

    double sxx = 0.0;
    double sxy = 0.0;
    double sxz = 0.0;
    double syx = 0.0;
    double syy = 0.0;
    double syz = 0.0;
    double szx = 0.0;
    double szy = 0.0;
    double szz = 0.0;

    for (std::size_t i = 0; i < p.size(); ++i) {
        sxx += p[i].x * q[i].x;
        sxy += p[i].x * q[i].y;
        sxz += p[i].x * q[i].z;
        syx += p[i].y * q[i].x;
        syy += p[i].y * q[i].y;
        syz += p[i].y * q[i].z;
        szx += p[i].z * q[i].x;
        szy += p[i].z * q[i].y;
        szz += p[i].z * q[i].z;
    }

    const double trace = sxx + syy + szz;
    const std::array<std::array<double, 4>, 4> horn{{
        {{trace, syz - szy, szx - sxz, sxy - syx}},
        {{syz - szy, sxx - syy - szz, sxy + syx, szx + sxz}},
        {{szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy}},
        {{sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz}},
    }};

    const double max_eval = largest_eigenvalue_symmetric_4x4(horn);
    const double squared = std::max(0.0, (norm_squared_sum(p) + norm_squared_sum(q) - 2.0 * max_eval) /
                                           static_cast<double>(p.size()));
    return std::sqrt(squared);
}

static Molecule reordered_molecule(const Molecule& target, const std::vector<std::size_t>& mapping) {
    Molecule out;
    out.symbols.reserve(mapping.size());
    out.coords.reserve(mapping.size());
    for (std::size_t target_index : mapping) {
        out.symbols.push_back(target.symbols[target_index]);
        out.coords.push_back(target.coords[target_index]);
    }
    return out;
}

static void search_reordered_rmsd(
    const Molecule& reference,
    const Molecule& target,
    const std::vector<std::vector<std::size_t>>& candidates,
    const std::vector<std::vector<bool>>& reference_bonded,
    const std::vector<std::vector<bool>>& target_bonded,
    std::size_t atom_index,
    std::vector<bool>& used,
    std::vector<std::size_t>& mapping,
    std::size_t max_mappings,
    RmsdMatch& match
) {
    if (match.mappings_checked >= max_mappings) {
        return;
    }

    if (atom_index == reference.coords.size()) {
        const Molecule reordered = reordered_molecule(target, mapping);
        const double rmsd = aligned_rmsd(reference, reordered);
        ++match.mappings_checked;
        if (rmsd < match.best) {
            match.second = match.best;
            match.best = rmsd;
            match.best_mapping = mapping;
        } else if (rmsd < match.second) {
            match.second = rmsd;
        }
        return;
    }

    for (std::size_t target_index : candidates[atom_index]) {
        if (used[target_index]) {
            continue;
        }

        // Preserve graph structure as the mapping is built. This prevents the
        // old too-broad behavior where any atom of the same element could swap.
        bool graph_consistent = true;
        for (std::size_t previous = 0; previous < atom_index; ++previous) {
            if (reference_bonded[atom_index][previous] != target_bonded[target_index][mapping[previous]]) {
                graph_consistent = false;
                break;
            }
        }
        if (!graph_consistent) {
            continue;
        }

        used[target_index] = true;
        mapping[atom_index] = target_index;
        search_reordered_rmsd(
            reference,
            target,
            candidates,
            reference_bonded,
            target_bonded,
            atom_index + 1,
            used,
            mapping,
            max_mappings,
            match
        );
        used[target_index] = false;
    }
}

static RmsdMatch best_graph_preserving_reordered_rmsd(
    const Molecule& reference,
    const Molecule& target,
    std::size_t max_mappings,
    double bond_scale
) {
    if (reference.coords.size() != target.coords.size()) {
        throw std::runtime_error("Atom count differs during reordered RMSD comparison");
    }

    // Reorder mode is graph-preserving, not merely element-preserving.
    // Candidate atoms must have matching local signatures and pass the
    // bond/non-bond consistency checks in search_reordered_rmsd.
    const MolecularGraph reference_graph = infer_graph(reference, bond_scale);
    const MolecularGraph target_graph = infer_graph(target, bond_scale);

    std::vector<std::vector<std::size_t>> candidates(reference.coords.size());
    for (std::size_t i = 0; i < reference.symbols.size(); ++i) {
        for (std::size_t j = 0; j < target.symbols.size(); ++j) {
            if (reference_graph.signatures[i] == target_graph.signatures[j]) {
                candidates[i].push_back(j);
            }
        }
        if (candidates[i].empty()) {
            return {};
        }
    }

    std::vector<std::size_t> order(reference.coords.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        return candidates[lhs].size() < candidates[rhs].size();
    });

    std::vector<std::vector<std::size_t>> ordered_candidates(reference.coords.size());
    std::vector<std::vector<bool>> ordered_reference_bonded(reference.coords.size(), std::vector<bool>(reference.coords.size(), false));
    Molecule ordered_reference;
    ordered_reference.symbols.reserve(reference.coords.size());
    ordered_reference.coords.reserve(reference.coords.size());
    for (std::size_t pos = 0; pos < order.size(); ++pos) {
        ordered_candidates[pos] = candidates[order[pos]];
        ordered_reference.symbols.push_back(reference.symbols[order[pos]]);
        ordered_reference.coords.push_back(reference.coords[order[pos]]);
    }
    for (std::size_t i = 0; i < order.size(); ++i) {
        for (std::size_t j = 0; j < order.size(); ++j) {
            ordered_reference_bonded[i][j] = reference_graph.bonded[order[i]][order[j]];
        }
    }

    std::vector<bool> used(target.coords.size(), false);
    std::vector<std::size_t> mapping(reference.coords.size());
    RmsdMatch match;
    search_reordered_rmsd(
        ordered_reference,
        target,
        ordered_candidates,
        ordered_reference_bonded,
        target_graph.bonded,
        0,
        used,
        mapping,
        max_mappings,
        match
    );
    if (!match.best_mapping.empty()) {
        std::vector<std::size_t> original_order_mapping(reference.coords.size());
        for (std::size_t pos = 0; pos < order.size(); ++pos) {
            original_order_mapping[order[pos]] = match.best_mapping[pos];
        }
        match.best_mapping = std::move(original_order_mapping);
    }
    return match;
}

static void assert_same_atom_indexing(const Molecule& reference, const Molecule& mol, const std::string& path) {
    if (reference.coords.size() != mol.coords.size()) {
        throw std::runtime_error("Atom count differs from first conformer in " + path);
    }
    for (std::size_t i = 0; i < reference.symbols.size(); ++i) {
        if (reference.symbols[i] != mol.symbols[i]) {
            std::ostringstream msg;
            msg << "Atom index mismatch in " << path << " at zero-based index " << i << ": expected "
                << reference.symbols[i] << " but found " << mol.symbols[i]
                << ". Reorder atoms before deduplication.";
            throw std::runtime_error(msg.str());
        }
    }
}

static void assert_same_formula(const Molecule& reference, const Molecule& mol, const std::string& path) {
    if (reference.coords.size() != mol.coords.size()) {
        throw std::runtime_error("Atom count differs from first conformer in " + path);
    }
    std::vector<std::string> a = reference.symbols;
    std::vector<std::string> b = mol.symbols;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    if (a != b) {
        throw std::runtime_error("Element counts differ from first conformer in " + path);
    }
}

class Conformer_Batch {
public:
    struct Remove_Duplicates_Result {
        std::vector<UniqueEntry> unique;
        std::vector<DuplicateEntry> duplicates;
    };

    static Conformer_Batch from_xyz_files(const std::vector<std::string>& paths) {
        // Loader 1a: explicit single-XYZ files. Each file contributes one
        // conformer record and keeps its path as the source label.
        Conformer_Batch batch;
        for (const std::string& path : paths) {
            std::ifstream in(path);
            if (!in) {
                throw std::runtime_error("Could not open XYZ file: " + path);
            }
            Molecule mol;
            std::string comment;
            if (!read_xyz_record(in, mol, comment)) {
                throw std::runtime_error("Empty XYZ file: " + path);
            }
            batch.records_.push_back({batch.records_.size(), path, comment, std::move(mol)});
        }
        return batch;
    }

    static Conformer_Batch from_xyz_directory(const std::filesystem::path& directory) {
        // Loader 1b: all .xyz files under a directory, sorted for reproducible
        // representative selection.
        if (!std::filesystem::is_directory(directory)) {
            throw std::runtime_error("Not a directory: " + directory.string());
        }

        std::vector<std::string> paths;
        for (const auto& entry : std::filesystem::directory_iterator(directory)) {
            if (!entry.is_regular_file()) {
                continue;
            }
            if (entry.path().extension() == ".xyz") {
                paths.push_back(entry.path().string());
            }
        }
        std::sort(paths.begin(), paths.end());
        if (paths.empty()) {
            throw std::runtime_error("No .xyz files found under directory: " + directory.string());
        }
        return from_xyz_files(paths);
    }

    static Conformer_Batch from_multi_xyz(const std::filesystem::path& path) {
        // Loader 2: concatenated multi-XYZ. Source labels use file#record_index
        // so duplicate reports can point back to the original record.
        std::ifstream in(path);
        if (!in) {
            throw std::runtime_error("Could not open multi-XYZ file: " + path.string());
        }

        Conformer_Batch batch;
        while (true) {
            Molecule mol;
            std::string comment;
            if (!read_xyz_record(in, mol, comment)) {
                break;
            }
            std::ostringstream source;
            source << path.string() << "#" << batch.records_.size();
            batch.records_.push_back({batch.records_.size(), source.str(), comment, std::move(mol)});
        }
        if (batch.records_.empty()) {
            throw std::runtime_error("No XYZ records found in multi-XYZ file: " + path.string());
        }
        return batch;
    }

    std::size_t size() const {
        return records_.size();
    }

    const std::vector<Conformer_Record>& records() const {
        return records_;
    }

    void index_cleanup(std::size_t max_mappings, double bond_scale, double ambiguity_gap) {
        // Reorder every conformer to match the first conformer's atom indexing.
        // This method is intentionally public because callers may want a cleaned
        // conformer set even when they are not removing duplicates yet.
        if (records_.empty()) {
            return;
        }

        const Molecule reference = records_.front().mol;
        for (std::size_t i = 0; i < records_.size(); ++i) {
            Conformer_Record& record = records_[i];
            assert_same_formula(reference, record.mol, record.source);
            if (i == 0) {
                continue;
            }

            const RmsdMatch match =
                best_graph_preserving_reordered_rmsd(reference, record.mol, max_mappings, bond_scale);
            validate_mapping_result(match, record.source, max_mappings, ambiguity_gap);
            record.mol = reordered_molecule(record.mol, match.best_mapping);
            record.comment = "index-cleaned from " + record.source;
        }
    }

    Remove_Duplicates_Result remove_duplicates(
        double tolerance,
        bool run_index_cleanup,
        std::size_t max_mappings,
        double bond_scale,
        double ambiguity_gap
    ) {
        // run_index_cleanup is the class-level option that lets duplicate
        // removal handle reordered XYZ files in one call.
        if (run_index_cleanup) {
            index_cleanup(max_mappings, bond_scale, ambiguity_gap);
        }

        Remove_Duplicates_Result result;
        if (records_.empty()) {
            return result;
        }

        const Molecule reference = records_.front().mol;
        for (const Conformer_Record& record : records_) {
            assert_same_atom_indexing(reference, record.mol, record.source);

            double best_rmsd = std::numeric_limits<double>::infinity();
            const UniqueEntry* best = nullptr;
            for (const UniqueEntry& candidate : result.unique) {
                const double rmsd = aligned_rmsd(candidate.mol, record.mol);
                if (rmsd < best_rmsd) {
                    best_rmsd = rmsd;
                    best = &candidate;
                }
            }

            if (best != nullptr && best_rmsd <= tolerance) {
                result.duplicates.push_back({record.input_index, record.source, best->input_index, best->path, best_rmsd});
            } else {
                result.unique.push_back({record.input_index, record.source, record.mol});
            }
        }
        return result;
    }

    void write_records(const std::filesystem::path& outdir, const std::string& prefix) const {
        std::filesystem::create_directories(outdir);
        for (std::size_t i = 0; i < records_.size(); ++i) {
            std::ostringstream name;
            name << prefix << "_" << std::setw(4) << std::setfill('0') << (i + 1) << ".xyz";
            write_xyz(outdir / name.str(), records_[i].mol, records_[i].comment);
        }
    }

private:
    static void validate_mapping_result(
        const RmsdMatch& match,
        const std::string& source,
        std::size_t max_mappings,
        double ambiguity_gap
    ) {
        if (match.mappings_checked >= max_mappings && !std::isfinite(match.best)) {
            throw std::runtime_error("Mapping limit reached before a valid mapping was checked for " + source);
        }
        if (match.mappings_checked >= max_mappings) {
            throw std::runtime_error("Mapping limit reached for " + source + "; increase --max-mappings or use RDKit graph matching");
        }
        if (match.best_mapping.empty()) {
            throw std::runtime_error("No graph-compatible atom mapping found for " + source);
        }
        (void)ambiguity_gap;
        // Highly symmetric molecules can have several equally good graph-valid
        // mappings. That is acceptable here: keep the lowest-RMSD mapping found.
    }

    std::vector<Conformer_Record> records_;
};

#ifndef CONFORMER_TOOLKIT_NO_MAIN
static void usage(const char* program) {
    std::cerr << "Usage: " << program
              << " [--tolerance angstrom] [--allow-reorder] [--max-mappings n] [--ambiguity-gap angstrom]"
              << " [--bond-scale scale] [--xyz-dir dir | --multi-xyz file] [--index-cleanup-only]"
              << " [--write-cleaned dir] [--write-unique dir] conformer1.xyz conformer2.xyz [...]\n"
              << "Keeps the first representative of each unique conformer group.\n"
              << "--xyz-dir loads all .xyz files under a directory, sorted by path.\n"
              << "--multi-xyz loads multiple XYZ records from one concatenated XYZ file.\n"
              << "By default all conformers must have the same atom symbols in the same order.\n"
              << "--allow-reorder runs Conformer_Batch::index_cleanup before duplicate removal.\n";
}

int main(int argc, char** argv) {
    double tolerance = 1e-3;
    double ambiguity_gap = 1e-6;
    double bond_scale = 1.1;
    bool allow_reorder = false;
    bool index_cleanup_only = false;
    std::size_t max_mappings = 1000000;
    std::filesystem::path xyz_dir;
    std::filesystem::path multi_xyz;
    std::filesystem::path write_cleaned_dir;
    std::filesystem::path write_unique_dir;
    std::vector<std::string> paths;

    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--tolerance") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                tolerance = std::stod(argv[++i]);
            } else if (arg == "--allow-reorder") {
                allow_reorder = true;
            } else if (arg == "--index-cleanup-only") {
                index_cleanup_only = true;
            } else if (arg == "--bond-scale") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                bond_scale = std::stod(argv[++i]);
            } else if (arg == "--xyz-dir") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                xyz_dir = argv[++i];
            } else if (arg == "--multi-xyz") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                multi_xyz = argv[++i];
            } else if (arg == "--max-mappings") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                max_mappings = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--ambiguity-gap") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                ambiguity_gap = std::stod(argv[++i]);
            } else if (arg == "--write-unique") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                write_unique_dir = argv[++i];
            } else if (arg == "--write-cleaned") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                write_cleaned_dir = argv[++i];
            } else if (arg == "-h" || arg == "--help") {
                usage(argv[0]);
                return 0;
            } else if (arg.rfind("--", 0) == 0) {
                usage(argv[0]);
                return 2;
            } else {
                paths.push_back(arg);
            }
        }

        const int input_mode_count = (!paths.empty() ? 1 : 0) + (!xyz_dir.empty() ? 1 : 0) + (!multi_xyz.empty() ? 1 : 0);
        if (input_mode_count != 1) {
            usage(argv[0]);
            std::cerr << "error: provide exactly one input mode: positional XYZ files, --xyz-dir, or --multi-xyz\n";
            return 2;
        }
        if (tolerance < 0.0) {
            throw std::runtime_error("--tolerance must be non-negative");
        }
        if (max_mappings == 0) {
            throw std::runtime_error("--max-mappings must be greater than zero");
        }
        if (bond_scale <= 0.0) {
            throw std::runtime_error("--bond-scale must be greater than zero");
        }

        Conformer_Batch batch = [&]() {
            if (!xyz_dir.empty()) {
                return Conformer_Batch::from_xyz_directory(xyz_dir);
            }
            if (!multi_xyz.empty()) {
                return Conformer_Batch::from_multi_xyz(multi_xyz);
            }
            return Conformer_Batch::from_xyz_files(paths);
        }();

        const auto start = std::chrono::steady_clock::now();
        if (index_cleanup_only) {
            batch.index_cleanup(max_mappings, bond_scale, ambiguity_gap);
            if (write_cleaned_dir.empty()) {
                throw std::runtime_error("--index-cleanup-only requires --write-cleaned");
            }
            batch.write_records(write_cleaned_dir, "cleaned");
        }
        const Conformer_Batch::Remove_Duplicates_Result result =
            index_cleanup_only
                ? Conformer_Batch::Remove_Duplicates_Result{}
                : batch.remove_duplicates(tolerance, allow_reorder, max_mappings, bond_scale, ambiguity_gap);
        const auto stop = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = stop - start;

        if (!write_cleaned_dir.empty() && !index_cleanup_only && allow_reorder) {
            batch.write_records(write_cleaned_dir, "cleaned");
        }

        if (!write_unique_dir.empty() && !index_cleanup_only) {
            std::filesystem::create_directories(write_unique_dir);
            for (std::size_t i = 0; i < result.unique.size(); ++i) {
                std::ostringstream name;
                name << "unique_" << std::setw(4) << std::setfill('0') << (i + 1) << ".xyz";
                write_xyz(write_unique_dir / name.str(), result.unique[i].mol, "representative from " + result.unique[i].path);
            }
        }

        std::cout << std::fixed << std::setprecision(10);
        std::cout << "input_count " << batch.size() << "\n";
        std::cout << "unique_count " << result.unique.size() << "\n";
        std::cout << "duplicate_count " << result.duplicates.size() << "\n";
        std::cout << "tolerance_angstrom " << tolerance << "\n";
        std::cout << "elapsed_seconds " << elapsed.count() << "\n";
        if (index_cleanup_only) {
            std::cout << "mode index_cleanup_only\n";
        }
        for (std::size_t i = 0; i < result.unique.size(); ++i) {
            std::cout << "unique " << i << " input_index " << result.unique[i].input_index << " path " << result.unique[i].path
                      << "\n";
        }
        for (const DuplicateEntry& duplicate : result.duplicates) {
            std::cout << "duplicate input_index " << duplicate.input_index << " path " << duplicate.path
                      << " representative_input_index " << duplicate.representative_input_index
                      << " representative_path " << duplicate.representative_path << " rmsd_angstrom "
                      << duplicate.rmsd << "\n";
        }

        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "error: " << exc.what() << "\n";
        return 2;
    }
}
#endif
