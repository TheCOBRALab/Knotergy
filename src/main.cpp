#include <iostream>
#include <string>
#include <unordered_set>

#include "energy/ComputeEnergy.hpp"
#include "io/OutputManager.hpp"
#include "io/PseudoknotParams.hpp"
#include "io/RNAInputManager.hpp"
#include "loop_tree/LoopFactory.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
#include "preprocessing/RNAEntry.hpp"

#define KNOTERGY_VERSION "0.1.3"

namespace {
void help() {
    std::cout
        << "Usage: ./Knotergy [options]\n"
        << "Options:\n"
        << "  -h, --help                            Show this help message\n"
        << "  -V, --version                         Print version and exit\n"
        << "  -v, --verbose                         Enable verbose output\n"
        << "  -s, --sequence <string>               Input sequence\n"
        << "  -r, --structure <string>              Input structure\n"
        << "  -i, --input <file>                    Input file\n"
        << "  -p, --paramFile <file>                Parameter file\n"
        << "  -k, --pk-paramFile <file>             Pseudoknot parameter file\n"
        << "  -m, --mod-params <none|path|file>     Directory containing modified base parameter files\n"
        << "  -e, --round                           Rounds all decimal places in pseudoknot calculations\n"
        << "  -d, --dangle                          Specify the dangle model to be used (base is 2)\n";
}

// cleans white space from arg
std::string get_trimmed_arg(int& i, int argc, char** argv) {
    if (i + 1 >= argc) return "";
    std::string value = argv[++i];
    knotergy::trim(value);
    return value;
}

// converts arg to int, with default value if conversion fails
int get_numerical_arg(int& i, int argc, char** argv, int default_value = 0) {
    std::string value_str = get_trimmed_arg(i, argc, argv);
    if (value_str.empty()) return default_value;
    try {
        return std::stoi(value_str);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid numerical argument: " << value_str << std::endl;
        return default_value;
    } catch (const std::out_of_range& e) {
        std::cerr << "Numerical argument out of range: " << value_str << std::endl;
        return default_value;
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::string sequence = "";
    std::string structure = "";
    std::string input_file = "";
    std::string output_file = "";
    std::string parameter_file = "";
    std::string modifications = "7I6P9D";
    std::vector<std::string> mod_param_paths; 
    const std::string default_mod_path = std::string(KNOTERGY_SOURCE_DIR) + "/params/modified_bases";
    std::string pseudo_param_file = std::string(KNOTERGY_SOURCE_DIR) + "/params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json";
    bool round = false;
    bool verbose = false;
    int dangle = 2;

    // ------------------------- Parse Through Flags -----------------------
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-s" || arg == "--sequence") && argc >= i + 1) {
            sequence = get_trimmed_arg(i, argc, argv);
        } else if ((arg == "-r" || arg == "--structure") && argc >= i + 1) {
            structure = get_trimmed_arg(i, argc, argv);
        } else if ((arg == "-i" || arg == "--input") && argc >= i + 1) {
            input_file = get_trimmed_arg(i, argc, argv);
        } else if ((arg == "-p" || arg == "--paramFile") && argc >= i + 1) {
            parameter_file = get_trimmed_arg(i, argc, argv);
        } else if (arg == "-e" || arg == "--round") {
            round = true;
        } else if (arg == "-m" || arg == "--mod-params") {
            // -m /path/to/params
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                mod_param_paths.push_back(get_trimmed_arg(i, argc, argv));
            } else { // -m with no path, use default
                mod_param_paths.push_back(default_mod_path);
            }
        } else if (arg == "-d" || arg == "--dangle") {
            dangle = get_numerical_arg(i, argc, argv, dangle);
        } else if ((arg == "-k" || arg == "--pk-paramFile") && argc >= i + 1) {
            pseudo_param_file = get_trimmed_arg(i, argc, argv);
        } else if (arg == "-h" || arg == "--help") {
            help();
            return 0;
        } else if (arg == "-V" || arg == "--version") {
            std::cout << "Knotergy " << KNOTERGY_VERSION << std::endl;
            return 0;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else {
            std::cerr << "Unknown option or missing value: " << arg << std::endl;
            return 1;
        }
    }

    // ------------------------- Validate Inputs -----------------------
    // Get sequence if not provided
    if (sequence.empty() && input_file.empty()) {
        std::cout << "Sequence : ";
        std::cin >> sequence;
        knotergy::trim(sequence);
    }

    // Get structure if not provided
    if (structure.empty() && input_file.empty()) {
        std::cout << "Structure: ";
        std::cin >> structure;
        knotergy::trim(structure);
    }

    knotergy::trim(input_file);
    if (!input_file.empty() && !knotergy::FileUtils::file_exists(input_file)) {
        std::cerr << "Input file not found: " << input_file << std::endl;
        return 1;
    }

    knotergy::trim(parameter_file);
    if (!parameter_file.empty() && !knotergy::FileUtils::file_exists(parameter_file)) {
        std::cerr << "Parameter file not found: " << parameter_file << std::endl;
        return 1;
    }

    for (std::string& mod_path : mod_param_paths) {
        knotergy::trim(mod_path);
        if (!knotergy::FileUtils::file_exists(mod_path)) {
            std::cerr << "Modified bases parameter path not found: " << mod_path << std::endl;
            return 1;
        }
    }

    if (structure.length() >= 2147483647) { // Prevent overflow of int in energy calculations
        std::cerr << "Error: Structure length exceeds maximum allowed size of 2,147,483,647" << std::endl;
        return 1;
    }

    // ------------------------- Load ViennaRNA Parameters -----------------------
    knotergy::vrna_md_param vp =
        knotergy::ViennaParams::load_energy_parameters(parameter_file, dangle, sequence);

    // ------------------------- Load Pseudoknot Parameters -----------------------
    knotergy::pk_param pkp = knotergy::PseudoknotParams::load_pk_param(pseudo_param_file);

    // ------------------------- Load Modified Base Parameters -----------------------
    std::vector<knotergy::modified_base_param> mp;
    for (const std::string& mod_path : mod_param_paths) {
        std::vector<knotergy::modified_base_param> additional_mp =
            knotergy::ViennaParams::load_modified_energy_parameters(mod_path);
        mp.insert(mp.end(), additional_mp.begin(), additional_mp.end());
    }

    //------------------------- Reading Inputs From File -----------------------------
    std::vector<knotergy::RNAEntry> inputs =
        knotergy::RNAInputManager::get_all_inputs(input_file, sequence, structure);

    // ------------------------- Print Parameter Report -----------------------
    knotergy::OutputManager::print_parameter_report(vp, pkp, mp);

    //------------------------- Main Processing Loop ----------------------------
    for (const knotergy::RNAEntry& rna : inputs) {
        std::cout << "\n--------- Name: " << rna.name << " ---------" << std::endl;

        // Preprocess the RNA entry to compute pairings, closed regions, etc.
        const knotergy::ProcessedRNAEntry& processed_rna = knotergy::RNAProcessor::process_rna(rna, mp);

        // Builds loop tree
        knotergy::LoopFactory factory(processed_rna);

        // Computes the energy
        knotergy::ComputeEnergy energy_calculator(factory.get_root_node(), processed_rna, vp, pkp,
                                                  mp, round, verbose);

        // Output results
        if (energy_calculator.getInfiniteEnergyFlag()) {
            printf("\nENERGY: Infinite (%.4f kcal/mol)\n", energy_calculator.getEnergy());
        } else {
            printf("\nENERGY: %.4f kcal/mol\n", energy_calculator.getEnergy());
        }
    }
    return 0;
}