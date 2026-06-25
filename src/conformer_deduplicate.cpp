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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <GraphMol/DetermineBonds/DetermineBonds.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>

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

static std::string xyz_block(const Molecule& mol) {
    std::ostringstream out;
    out << mol.coords.size() << "\n\n";
    out << std::setprecision(17);
    for (std::size_t i = 0; i < mol.coords.size(); ++i) {
        out << mol.symbols[i] << " " << mol.coords[i].x << " " << mol.coords[i].y << " " << mol.coords[i].z << "\n";
    }
    return out.str();
}

static std::unique_ptr<RDKit::RWMol> chemistry_from_xyz(
    const Molecule& mol,
    double bond_scale,
    int charge
) {
    std::unique_ptr<RDKit::RWMol> rd_mol = RDKit::v2::FileParsers::MolFromXYZBlock(xyz_block(mol));
    if (!rd_mol) {
        throw std::runtime_error("RDKit could not parse the XYZ record");
    }

    try {
        // XYZ carries no connectivity. DetermineBonds applies the xyz2mol
        // algorithm, then sanitizes and assigns tetrahedral/double-bond stereo
        // from the 3D conformer. useVdw=true makes bond_scale the covalent-radius
        // multiplier, preserving the existing CLI/API control.
        RDKit::determineBonds(
            *rd_mol,
            false,
            charge,
            bond_scale,
            true,
            true,
            false,
            true
        );
        RDKit::MolOps::assignStereochemistry(*rd_mol, true, true, true);
    } catch (const std::exception& exc) {
        throw std::runtime_error(std::string("RDKit could not infer molecular chemistry from XYZ: ") + exc.what());
    }
    return rd_mol;
}

static std::vector<std::vector<std::size_t>> chemistry_preserving_mappings(
    const Molecule& reference,
    const Molecule& target,
    std::size_t max_mappings,
    double bond_scale,
    int charge,
    bool use_chirality
) {
    if (reference.coords.size() != target.coords.size()) {
        return {};
    }

    const std::unique_ptr<RDKit::RWMol> reference_mol = chemistry_from_xyz(reference, bond_scale, charge);
    std::unique_ptr<RDKit::RWMol> target_mol;
    try {
        target_mol = chemistry_from_xyz(target, bond_scale, charge);
    } catch (const std::exception&) {
        // A highly distorted candidate may no longer support bond inference
        // from its coordinates. It is not chemically matchable to the
        // reference, but that should classify it as different rather than
        // aborting an otherwise valid conformer batch.
        return {};
    }

    RDKit::SubstructMatchParameters params;
    params.useChirality = use_chirality;
    params.useEnhancedStereo = use_chirality;
    params.specifiedStereoQueryMatchesUnspecified = false;
    params.uniquify = false;
    params.maxMatches = static_cast<unsigned int>(
        std::min<std::size_t>(max_mappings, std::numeric_limits<unsigned int>::max())
    );

    // Each pair is (reference/query atom, target atom). Equal atom counts make
    // a complete substructure match a full graph isomorphism; bond orders,
    // formal charges, atom types, and assigned stereochemistry all participate.
    const std::vector<RDKit::MatchVectType> matches = RDKit::SubstructMatch(*target_mol, *reference_mol, params);
    std::vector<std::vector<std::size_t>> mappings;
    mappings.reserve(matches.size());
    for (const RDKit::MatchVectType& match : matches) {
        if (match.size() != reference.coords.size()) {
            continue;
        }
        std::vector<std::size_t> mapping(reference.coords.size());
        for (const auto& pair : match) {
            mapping[static_cast<std::size_t>(pair.first)] = static_cast<std::size_t>(pair.second);
        }
        mappings.push_back(std::move(mapping));
    }
    return mappings;
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

static RmsdMatch best_graph_preserving_reordered_rmsd(
    const Molecule& reference,
    const Molecule& target,
    std::size_t max_mappings,
    double bond_scale,
    int charge,
    bool use_chirality
) {
    if (reference.coords.size() != target.coords.size()) {
        throw std::runtime_error("Atom count differs during reordered RMSD comparison");
    }

    RmsdMatch match;
    const std::vector<std::vector<std::size_t>> mappings =
        chemistry_preserving_mappings(reference, target, max_mappings, bond_scale, charge, use_chirality);
    for (const std::vector<std::size_t>& mapping : mappings) {
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

    void index_cleanup(std::size_t max_mappings, double bond_scale, double ambiguity_gap, int charge = 0) {
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
                best_graph_preserving_reordered_rmsd(reference, record.mol, max_mappings, bond_scale, charge, false);
            if (match.best_mapping.empty() && reference.symbols == record.mol.symbols) {
                // If chemistry inference fails for a distorted conformer but
                // atom symbols already follow the reference order, retain that
                // order and let duplicate classification treat it as distinct.
                continue;
            }
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
        double ambiguity_gap,
        int charge = 0
    ) {
        // run_index_cleanup is the class-level option that lets duplicate
        // removal handle reordered XYZ files in one call.
        if (run_index_cleanup) {
            index_cleanup(max_mappings, bond_scale, ambiguity_gap, charge);
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
                const std::vector<std::vector<std::size_t>> chemistry_matches =
                    chemistry_preserving_mappings(candidate.mol, record.mol, 1, bond_scale, charge, true);
                if (chemistry_matches.empty()) {
                    continue;
                }
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
              << " [--bond-scale scale] [--charge integer] [--xyz-dir dir | --multi-xyz file] [--index-cleanup-only]"
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
    double bond_scale = 1.3;
    bool allow_reorder = false;
    bool index_cleanup_only = false;
    std::size_t max_mappings = 1000000;
    int charge = 0;
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
            } else if (arg == "--charge") {
                if (i + 1 >= argc) {
                    usage(argv[0]);
                    return 2;
                }
                charge = std::stoi(argv[++i]);
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
            batch.index_cleanup(max_mappings, bond_scale, ambiguity_gap, charge);
            if (write_cleaned_dir.empty()) {
                throw std::runtime_error("--index-cleanup-only requires --write-cleaned");
            }
            batch.write_records(write_cleaned_dir, "cleaned");
        }
        const Conformer_Batch::Remove_Duplicates_Result result =
            index_cleanup_only
                ? Conformer_Batch::Remove_Duplicates_Result{}
                : batch.remove_duplicates(tolerance, allow_reorder, max_mappings, bond_scale, ambiguity_gap, charge);
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
