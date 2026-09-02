#include "io/input/CliArgs.hpp"

#include "io/output/OutputManager.hpp"
#include "io/output/colors.hpp"
#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "utils/FileUtils.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace knotergy {

namespace {

// Largest structure length we can index with a signed 32-bit int.
constexpr std::size_t MAX_STRUCTURE_LENGTH =
    static_cast<std::size_t>(std::numeric_limits<int>::max());

// True when next argument is a value rather than another flag (e.g. "-s" or "--sequence").
bool next_is_value(int i, int argc, char** argv) { return i + 1 < argc && argv[i + 1][0] != '-'; }

// Reads the value following argv[i], advancing i past it and trimming
// whitespace. Reports and returns false when the option has no value.
bool take_value(int& i, int argc, char** argv, std::string& out) {
    if (i + 1 >= argc) {
        std::cerr << ERROR << " Missing value for option: " << argv[i] << '\n';
        return false;
    }

    out = argv[++i];
    trim(out);
    return true;
}

// As take_value, but converts the value to an int.
bool take_int_value(int& i, int argc, char** argv, int& out) {
    std::string value;
    if (!take_value(i, argc, argv, value)) return false;

    try {
        out = std::stoi(value);
        return true;
    } catch (const std::invalid_argument&) {
        std::cerr << ERROR << " Invalid numerical argument: " << value << '\n';
    } catch (const std::out_of_range&) {
        std::cerr << ERROR << " Numerical argument out of range: " << value << '\n';
    }

    return false;
}

// As take_value, but converts the value to a double. `what` names the option
// in the error message, e.g. "temperature".
bool take_double_value(int& i, int argc, char** argv, const char* what, double& out) {
    std::string value;
    if (!take_value(i, argc, argv, value)) return false;

    try {
        out = std::stod(value);
        return true;
    } catch (const std::exception&) {
        std::cerr << ERROR << " Invalid " << what << " value: " << value << '\n';
        return false;
    }
}

// Parses digits glued onto a short flag, e.g. the "2" in "-d2".
bool parse_attached_int(const std::string& arg, const char* what, int& out) {
    try {
        out = std::stoi(arg.substr(2));
        return true;
    } catch (const std::exception&) {
        std::cerr << ERROR << " Invalid " << what << " value: " << arg.substr(2) << '\n';
        return false;
    }
}

// True when `arg` is a short flag with its value attached, e.g. "-d2".
bool is_attached_form(const std::string& arg, const char* flag) {
    return arg.size() > 2 && arg.rfind(flag, 0) == 0;
}

}  // namespace

// ------------------------------- Help -----------------------------------

void print_help() {
    std::cout << "Usage: ./Knotergy [options]\n"
              << "Options:\n"
              << "  -h, --help                            Show this help message\n"
              << "  -V, --version                         Print version and exit\n"
              << "  -v, --verbose                         Enable verbose output\n"
              << "  -s, --sequence <string>               Input sequence\n"
              << "  -r, --structure <string>              Input structure\n"
              << "  -i, --input <file>                    Input file\n"
              << "  -P, --paramFile <file>                Parameter file\n"
              << "  -k, --pk-paramFile <file>             Pseudoknot parameter file\n"
              << "  -m, --mod-params <none|path|file>     Directory containing modified base "
                 "parameter files\n"
              << "  -T, --temp <temperature>              Temperature in Celsius (default: 37)\n"
              << "  -e, --round                           Rounds all decimal places in pseudoknot "
                 "calculations\n"
              << "  -d, --dangles                         Specify the dangle model to be used "
                 "(base is 2)\n"
              << "      --pk-dangles                      Enable pseudoknot dangle calculations\n"
              << "      --show-input                      Show the input sequence and structure\n"
              << "      --disable-cache                   Disable parameter caching\n"
              << "      --salt <salt>                     Salt Concentration (default: 1.021)"
        //   << "      --convertU                        Convert T to U"
        //   << "      --convertT                        Convert U to T"
        //   << "      --efn2-correction                 Apply efn2 single-bulge correction\n"
        ;
}

bool has_verbose_flag(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-v" || arg == "--verbose") {
            return true;
        }
    }

    return false;
}

// ------------------------------ Parsing ---------------------------------

ParseStatus CliArgs::parse(int argc, char** argv, CliArgs& out) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "-s" || arg == "--sequence") {
            if (!take_value(i, argc, argv, out.sequence)) return ParseStatus::ExitFailure;

        } else if (arg == "-r" || arg == "--structure") {
            if (!take_value(i, argc, argv, out.structure)) return ParseStatus::ExitFailure;

        } else if (arg == "-i" || arg == "--input") {
            if (!take_value(i, argc, argv, out.input_file)) return ParseStatus::ExitFailure;

        } else if (arg == "-P" || arg == "--paramFile") {
            if (!take_value(i, argc, argv, out.vienna_param_file)) return ParseStatus::ExitFailure;

        } else if (arg == "-k" || arg == "--pk-paramFile") {
            if (!take_value(i, argc, argv, out.pseudo_param_file)) return ParseStatus::ExitFailure;

        } else if (arg == "-e" || arg == "--round") {
            // -e <n>, or bare -e for the default rounding method.
            if (next_is_value(i, argc, argv)) {
                if (!take_int_value(i, argc, argv, out.round_value))
                    return ParseStatus::ExitFailure;
            } else {
                out.round_value = static_cast<int>(RoundMethod::Bankers);
            }

        } else if (is_attached_form(arg, "-e")) {
            // Supports: -e2 (no space between -e and the number)
            if (!parse_attached_int(arg, "round", out.round_value)) return ParseStatus::ExitFailure;

        } else if (arg == "-d" || arg == "--dangles") {
            if (!take_int_value(i, argc, argv, out.dangle)) return ParseStatus::ExitFailure;

        } else if (is_attached_form(arg, "-d")) {
            // Supports: -d2 (no space between -d and the number)
            if (!parse_attached_int(arg, "dangle", out.dangle)) return ParseStatus::ExitFailure;

        } else if (arg == "-m" || arg == "--mod-params") {
            // -m /path/to/params, or bare -m for the default location.
            if (next_is_value(i, argc, argv)) {
                std::string mod_path;
                if (!take_value(i, argc, argv, mod_path)) return ParseStatus::ExitFailure;
                out.mod_param_paths.push_back(mod_path);
            } else {
                out.mod_param_paths.push_back(default_mod_param_path());
            }

        } else if (arg == "-T" || arg == "--temp") {
            if (!take_double_value(i, argc, argv, "temperature", out.temperature)) {
                return ParseStatus::ExitFailure;
            }

        } else if (arg == "--salt") {
            if (!take_double_value(i, argc, argv, "salt", out.salt)) {
                return ParseStatus::ExitFailure;
            }

        } else if (arg == "-h" || arg == "--help") {
            print_help();
            return ParseStatus::ExitSuccess;

        } else if (arg == "-V" || arg == "--version") {
            OutputManager::print_banner();
            return ParseStatus::ExitSuccess;

        } else if (arg == "-v" || arg == "--verbose") {
            out.verbose = true;

        } else if (arg == "--efn2-correction") {
            out.efn2_correction = true;

        } else if (arg == "--show-input") {
            out.show_input = true;

        } else if (arg == "--pk-dangles") {
            out.pk_dangles = true;

        } else if (arg == "--disable-cache") {
            out.disable_cache = true;

        } else if (arg == "--convertU") {
            out.convert = Convert::UtoT;

        } else if (arg == "--convertT") {
            out.convert = Convert::TtoU;

        } else {
            std::cerr << ERROR << " Unknown option: " << arg << '\n';
            return ParseStatus::ExitFailure;
        }
    }

    return ParseStatus::Ok;
}

// ---------------------------- Interactive -------------------------------

void CliArgs::prompt_for_missing_input() {
    if (!input_file.empty()) return;

    if (sequence.empty()) {
        std::cout << "Sequence : ";
        std::cin >> sequence;
        trim(sequence);
    }

    if (structure.empty()) {
        std::cout << "Structure: ";
        std::cin >> structure;
        trim(structure);
    }
}

// ---------------------------- Validation --------------------------------

bool CliArgs::validate() const {
    if (!input_file.empty() && !FileUtils::file_exists(input_file)) {
        std::cerr << ERROR << " Input file not found: " << input_file << '\n';
        return false;
    }

    if (!vienna_param_file.empty() && !FileUtils::file_exists(vienna_param_file)) {
        std::cerr << ERROR << " Parameter file not found: " << vienna_param_file << '\n';
        return false;
    }

    for (const std::string& mod_path : mod_param_paths) {
        if (!FileUtils::file_exists(mod_path)) {
            std::cerr << ERROR << " Modified bases parameter path not found: " << mod_path << '\n';
            return false;
        }
    }

    if (structure.length() > MAX_STRUCTURE_LENGTH) {
        std::cerr << ERROR << " Structure length exceeds maximum allowed size of "
                  << MAX_STRUCTURE_LENGTH << "\n";
        return false;
    }

    if (dangle < 0 || dangle > 3) {
        std::cerr << ERROR << " Invalid dangle value: " << dangle
                  << ". Dangle must be between 0 and 3.\n";
        return false;
    }

    if (round_value < 0 || round_value > 5) {
        std::cerr << ERROR << " Invalid round value: " << round_value
                  << ". Round must be between 0 and 5.\n";
        return false;
    }

    return true;
}

}  // namespace knotergy
