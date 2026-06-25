#define CONFORMER_TOOLKIT_NO_MAIN
#include "conformer_deduplicate.cpp"

static void identical_usage(const char* program) {
    std::cerr << "Usage: " << program
              << " conformer_a.xyz conformer_b.xyz [tolerance_angstrom]"
              << " [--repeat n] [--bond-scale scale] [--charge integer]\n"
              << "Default tolerance: 1e-3 Angstrom RMSD after optimal proper rotation.\n"
              << "RDKit infers and verifies connectivity, bond orders, charge, and stereochemistry.\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        identical_usage(argv[0]);
        return 2;
    }

    double tolerance = 1e-3;
    double bond_scale = 1.3;
    int charge = 0;
    std::size_t repeat = 1;

    try {
        for (int i = 3; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--repeat") {
                if (i + 1 >= argc) {
                    identical_usage(argv[0]);
                    return 2;
                }
                repeat = static_cast<std::size_t>(std::stoull(argv[++i]));
                if (repeat == 0) {
                    throw std::runtime_error("--repeat must be greater than zero");
                }
            } else if (arg == "--bond-scale") {
                if (i + 1 >= argc) {
                    identical_usage(argv[0]);
                    return 2;
                }
                bond_scale = std::stod(argv[++i]);
            } else if (arg == "--charge") {
                if (i + 1 >= argc) {
                    identical_usage(argv[0]);
                    return 2;
                }
                charge = std::stoi(argv[++i]);
            } else if (arg.rfind("--", 0) == 0) {
                identical_usage(argv[0]);
                return 2;
            } else {
                tolerance = std::stod(arg);
            }
        }
        if (tolerance < 0.0) {
            throw std::runtime_error("tolerance must be non-negative");
        }
        if (bond_scale <= 0.0) {
            throw std::runtime_error("--bond-scale must be greater than zero");
        }

        const Conformer_Batch a_batch = Conformer_Batch::from_xyz_files({argv[1]});
        const Conformer_Batch b_batch = Conformer_Batch::from_xyz_files({argv[2]});
        const Molecule& a = a_batch.records().front().mol;
        const Molecule& b = b_batch.records().front().mol;

        if (a.coords.size() != b.coords.size()) {
            std::cout << "different\nreason atom_count\n";
            return 1;
        }
        if (a.symbols != b.symbols) {
            std::cout << "different\nreason atom_symbols_or_order\n";
            return 1;
        }
        if (chemistry_preserving_mappings(a, b, 1, bond_scale, charge, true).empty()) {
            std::cout << "different\nreason chemistry_or_stereochemistry\n";
            return 1;
        }

        double rmsd = 0.0;
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < repeat; ++i) {
            rmsd = aligned_rmsd(a, b);
        }
        const auto stop = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = stop - start;

        std::cout << std::fixed << std::setprecision(10);
        std::cout << (rmsd <= tolerance ? "identical\n" : "different\n");
        std::cout << "rmsd_angstrom " << rmsd << "\n";
        std::cout << "tolerance_angstrom " << tolerance << "\n";
        if (repeat > 1) {
            std::cout << "repeat " << repeat << "\n";
            std::cout << "elapsed_seconds " << elapsed.count() << "\n";
            std::cout << "avg_microseconds " << (elapsed.count() * 1e6 / static_cast<double>(repeat)) << "\n";
        }
        return rmsd <= tolerance ? 0 : 1;
    } catch (const std::exception& exc) {
        std::cerr << "error: " << exc.what() << "\n";
        return 2;
    }
}
