#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
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

static Molecule read_xyz(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Could not open XYZ file: " + path);
    }

    std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("Empty XYZ file: " + path);
    }

    std::istringstream count_stream(line);
    std::size_t atom_count = 0;
    if (!(count_stream >> atom_count)) {
        throw std::runtime_error("First XYZ line is not an atom count: " + path);
    }

    std::getline(in, line);  // comment line

    Molecule mol;
    mol.symbols.reserve(atom_count);
    mol.coords.reserve(atom_count);

    for (std::size_t i = 0; i < atom_count; ++i) {
        if (!std::getline(in, line)) {
            throw std::runtime_error("XYZ ended before all atoms were read: " + path);
        }

        std::istringstream atom_stream(line);
        std::string symbol;
        Vec3 coord;
        if (!(atom_stream >> symbol >> coord.x >> coord.y >> coord.z)) {
            throw std::runtime_error("Bad XYZ atom line in " + path + ": " + line);
        }

        mol.symbols.push_back(symbol);
        mol.coords.push_back(coord);
    }

    return mol;
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

static void usage(const char* program) {
    std::cerr << "Usage: " << program << " conformer_a.xyz conformer_b.xyz [tolerance_angstrom] [--repeat n]\n"
              << "Default tolerance: 1e-3 Angstrom RMSD after optimal rigid alignment.\n"
              << "--repeat runs the comparison n times in-process and reports average time.\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 2;
    }

    double tolerance = 1e-3;
    std::size_t repeat = 1;

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--repeat") {
            if (i + 1 >= argc) {
                usage(argv[0]);
                return 2;
            }
            repeat = static_cast<std::size_t>(std::stoull(argv[++i]));
            if (repeat == 0) {
                throw std::runtime_error("--repeat must be greater than zero");
            }
        } else if (arg.rfind("--", 0) == 0) {
            usage(argv[0]);
            return 2;
        } else {
            tolerance = std::stod(arg);
        }
    }

    try {
        const Molecule a = read_xyz(argv[1]);
        const Molecule b = read_xyz(argv[2]);

        if (a.coords.size() != b.coords.size()) {
            std::cout << "different\nreason atom_count\n";
            return 1;
        }

        if (a.symbols != b.symbols) {
            std::cout << "different\nreason atom_symbols_or_order\n";
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
