#include "energy/ComputeEnergy.hpp"
#include "io/input/RNAInputManager.hpp"
#include "io/output/OutputManager.hpp"
#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopFactory.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
#include "preprocessing/RNAEntry.hpp"
#include "utils/colors.hpp"
#include "utils/common.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

void help() {
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
              << "  -e, --round                           Rounds all decimal places in pseudoknot "
                 "calculations\n"
              << "  -d, --dangles                         Specify the dangle model to be used "
                 "(base is 2)\n";
}

bool has_verbose_flag(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-v" || arg == "--verbose") {
            return true;
        }
    }

    return false;
}

// Cleans white space from arg.
std::string get_trimmed_arg(int& i, int argc, char** argv) {
    if (i + 1 >= argc) return "";

    std::string value = argv[++i];
    knotergy::trim(value);
    return value;
}

// Converts arg to int, with default value if conversion fails.
int get_numerical_arg(int& i, int argc, char** argv, int default_value = 0) {
    std::string value_str = get_trimmed_arg(i, argc, argv);

    if (value_str.empty()) return default_value;

    try {
        return std::stoi(value_str);
    } catch (const std::invalid_argument&) {
        std::cerr << ERROR << " Invalid numerical argument: " << value_str << '\n';
        return default_value;
    } catch (const std::out_of_range&) {
        std::cerr << ERROR << " Numerical argument out of range: " << value_str << '\n';
        return default_value;
    }
}

int run_knotergy(int argc, char** argv) {
    std::string sequence = "";
    std::string structure = "";
    std::string input_file = "";
    std::string output_file = "";
    std::string vienna_param_file = "";
    std::string pseudo_param_file = knotergy::default_pk_param_path();
    std::vector<std::string> mod_param_paths;
    std::string modifications = "7I6P9D";
    const bool use_color = knotergy::should_use_color();
    int round_value = 0;  // Default round value is 0 (no rounding)
    bool verbose = false;
    int dangle = 2;

    // ------------------------- Parse Through Flags -----------------------
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "-s" || arg == "--sequence") && i + 1 < argc) {
            sequence = get_trimmed_arg(i, argc, argv);

        } else if ((arg == "-r" || arg == "--structure") && i + 1 < argc) {
            structure = get_trimmed_arg(i, argc, argv);

        } else if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            input_file = get_trimmed_arg(i, argc, argv);

        } else if ((arg == "-P" || arg == "--paramFile") && i + 1 < argc) {
            vienna_param_file = get_trimmed_arg(i, argc, argv);

        } else if (arg == "-p") {
            std::cerr << WARNING << "-p is deprecated. " << "Use -P instead. "
                      << "-p will stop working on full release." << '\n';

            vienna_param_file = get_trimmed_arg(i, argc, argv);

        } else if (arg == "-e" || arg == "--round") {
            round_value =
                static_cast<int>(knotergy::RoundMethod::Bankers);  // Default to Banker's rounding
            // -m /path/to/params
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                round_value = get_numerical_arg(i, argc, argv, round_value);
            } else {
                // -e with no path, use default
                round_value = 1;
            }
        } else if (arg.rfind("-e", 0) == 0 && arg.size() > 2) {
            // Supports: -e2
            try {
                round_value = std::stoi(arg.substr(2));
            } catch (const std::exception&) {
                std::cerr << ERROR << " Invalid round value: " << arg.substr(2) << '\n';
                return EXIT_FAILURE;
            }
        } else if (arg == "-m" || arg == "--mod-params") {
            // -m /path/to/params
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                mod_param_paths.push_back(get_trimmed_arg(i, argc, argv));
            } else {
                // -m with no path, use default
                mod_param_paths.push_back(knotergy::default_mod_param_path());
            }

        } else if (arg == "-d" || arg == "--dangles") {
            dangle = get_numerical_arg(i, argc, argv, dangle);

        } else if (arg.rfind("-d", 0) == 0 && arg.size() > 2) {
            // Supports: -d2
            try {
                dangle = std::stoi(arg.substr(2));
            } catch (const std::exception&) {
                std::cerr << ERROR << " Invalid dangle value: " << arg.substr(2) << '\n';
                return EXIT_FAILURE;
            }

        } else if ((arg == "-k" || arg == "--pk-paramFile") && i + 1 < argc) {
            pseudo_param_file = get_trimmed_arg(i, argc, argv);

        } else if (arg == "-h" || arg == "--help") {
            help();
            return EXIT_SUCCESS;

        } else if (arg == "-V" || arg == "--version") {
            knotergy::OutputManager::print_banner();
            return EXIT_SUCCESS;

        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;

        } else {
            std::cerr << ERROR << " Unknown option or missing value: " << arg << '\n';
            return EXIT_FAILURE;
        }
    }

    // ------------------------- Validate Inputs -----------------------

    // Get sequence if not provided.
    if (sequence.empty() && input_file.empty()) {
        std::cout << "Sequence : ";
        std::cin >> sequence;
        knotergy::trim(sequence);
    }

    // Get structure if not provided.
    if (structure.empty() && input_file.empty()) {
        std::cout << "Structure: ";
        std::cin >> structure;
        knotergy::trim(structure);
    }

    knotergy::trim(input_file);
    if (!input_file.empty() && !knotergy::FileUtils::file_exists(input_file)) {
        std::cerr << ERROR << " Input file not found: " << input_file << '\n';
        return EXIT_FAILURE;
    }

    knotergy::trim(vienna_param_file);
    if (!vienna_param_file.empty() && !knotergy::FileUtils::file_exists(vienna_param_file)) {
        std::cerr << ERROR << " Parameter file not found: " << vienna_param_file << '\n';
        return EXIT_FAILURE;
    }

    for (std::string& mod_path : mod_param_paths) {
        knotergy::trim(mod_path);

        if (!knotergy::FileUtils::file_exists(mod_path)) {
            std::cerr << ERROR << " Modified bases parameter path not found: " << mod_path << '\n';
            return EXIT_FAILURE;
        }
    }

    if (structure.length() >= 2147483647) {
        std::cerr << ERROR << " Structure length exceeds maximum allowed size of 2,147,483,647\n";
        return EXIT_FAILURE;
    }

    if (dangle < 0 || dangle > 3) {
        std::cerr << ERROR << " Invalid dangle value: " << dangle
                  << ". Dangle must be 0, 1, 2, or 3.\n";
        return EXIT_FAILURE;
    }

    if (round_value < 0 || round_value > 5) {
        std::cerr << ERROR << " Invalid round value: " << round_value
                  << ". Round must be 0, 1, 2, 3, 4, or 5.\n";
        return EXIT_FAILURE;
    }
    const knotergy::RoundMethod round_method = static_cast<knotergy::RoundMethod>(round_value);

    // ------------------------- Load ViennaRNA Parameters -----------------------
    knotergy::vrna_md_param vp =
        knotergy::ViennaParams::load_energy_parameters(vienna_param_file, dangle, sequence);

    // ------------------------- Load Pseudoknot Parameters -----------------------
    knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pseudo_param_file, round_method);

    // ------------------------- Load Modified Base Parameters -----------------------
    std::vector<knotergy::modified_base_param> modified_params;

    for (const std::string& mod_path : mod_param_paths) {
        std::vector<knotergy::modified_base_param> additional_mp =
            knotergy::ModParams::load_modified_energy_parameters(mod_path);

        modified_params.insert(modified_params.end(), additional_mp.begin(), additional_mp.end());
    }

    knotergy::all_mod_params mp(modified_params);

    // ------------------------- Reading Inputs From File -----------------------------
    std::vector<knotergy::RNAEntry> inputs =
        knotergy::RNAInputManager::get_all_inputs(input_file, sequence, structure);

    // ------------------------- Print Parameter Report -----------------------
    knotergy::OutputManager::print_parameter_report(vp, pkp, mp);

    // ------------------------- Main Processing Loop ----------------------------
    for (const knotergy::RNAEntry& rna : inputs) {
        std::cout << "\n--------- Name: " << rna.name << " ---------" << '\n';

        // Preprocess the RNA entry to compute pair_table, closed regions, etc.
        const knotergy::ProcessedRNAEntry& processed_rna =
            knotergy::RNAProcessor::process_rna(rna, mp);

        // Build loop tree.
        knotergy::LoopFactory factory(processed_rna, vp);

        // Compute the energy.
        knotergy::ComputeEnergy energy_calculator(factory.get_root_node(), processed_rna, vp, pkp,
                                                  mp, verbose);

        // Output results.
        if (energy_calculator.getInfiniteEnergyFlag()) {
            printf("\nENERGY: Infinite (%.4f kcal/mol)\n", energy_calculator.getEnergy());
        } else if (use_color) {
            printf("\nENERGY:%s %.4f kcal/mol%s\n", ANSI_COLOR_GREEN, energy_calculator.getEnergy(),
                   ANSI_COLOR_RESET);
        } else {
            printf("\nENERGY: %.4f kcal/mol\n", energy_calculator.getEnergy());
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    const bool verbose = has_verbose_flag(argc, argv);

    try {
        return run_knotergy(argc, argv);

    } catch (const knotergy::DetailedException& e) {
        if (verbose) {
            std::cerr << e.detailed_message() << '\n';
        } else {
            std::cerr << ERROR << " " << e.what() << '\n';
        }

        return EXIT_FAILURE;

    } catch (const std::exception& e) {
        std::cerr << ERROR << " " << e.what() << '\n';
        return EXIT_FAILURE;

    } catch (...) {
        std::cerr << ERROR << " Unknown fatal error\n";
        return EXIT_FAILURE;
    }
}