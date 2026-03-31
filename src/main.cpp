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
        << "  -h, --help                    Show this help message\n"
        << "  -V, --version                 Print version and exit\n"
        << "  -v, --verbose                 Enable verbose output\n"
        << "  -s, --sequence <string>       RNA sequence\n"
        << "  -r, --structure <string>      Input structure\n"
        << "  -i, --input <file>            Input file\n"
        << "  -o, --output <file>           Output file\n"
        << "  -p, --paramFile <file>        Parameter file\n"
        << "  -k, --pk-paramFile <file>     Pseudoknot parameter file\n"
        << "  -m, --mod-dir <path|file>     Directory containing modified base parameter files\n"
        << "  -e, --round                   Rounds all decimal places in pseudoknot calculations\n"
        << "  -d, --dangle                  Specify the dangle model to be used (base is 2)\n";
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
    std::string mod_param_path = "./params/modified_bases";
    std::string pseudo_param_file = "./params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json";
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
        } else if ((arg == "-o" || arg == "--output") && argc >= i + 1) {
            output_file = get_trimmed_arg(i, argc, argv);
        } else if ((arg == "-p" || arg == "--paramFile") && argc >= i + 1) {
            parameter_file = get_trimmed_arg(i, argc, argv);
        } else if (arg == "-e" || arg == "--round") {
            round = true;
        } else if (arg == "-m" || arg == "--mod-file") {
            mod_param_path = get_trimmed_arg(i, argc, argv);
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

    knotergy::trim(mod_param_path);
    if (!mod_param_path.empty() && !knotergy::FileUtils::file_exists(mod_param_path)) {
        std::cerr << "Modified bases parameter file not found: " << mod_param_path << std::endl;
        return 1;
    }

    // ------------------------- Load ViennaRNA Parameters -----------------------
    knotergy::vrna_md_param vp =
        knotergy::ViennaParams::load_energy_parameters(parameter_file, dangle, sequence);

    // ------------------------- Load Pseudoknot Parameters -----------------------
    knotergy::pk_param pkp = knotergy::PseudoknotParams::load_pk_param(pseudo_param_file);

    // ------------------------- Load Modified Base Parameters -----------------------
    std::vector<knotergy::modified_base_param> mp =
        knotergy::ViennaParams::load_modified_energy_parameters(mod_param_path);

    //------------------------- Pre-processing and reading from files -----------------------------
    std::vector<knotergy::RNAEntry> inputs =
        knotergy::RNAInputManager::get_all_inputs(input_file, sequence, structure);
    std::vector<knotergy::ProcessedRNAEntry> processed_inputs =
        knotergy::RNAInputManager::process_inputs(inputs, mp);

    knotergy::OutputManager::print_parameter_report(vp, pkp, mp, verbose);

    //------------------------- Main Processing Loop ----------------------------
    for (const knotergy::ProcessedRNAEntry& processed_rna : processed_inputs) {
        std::cout << "\n--------- Name: " << processed_rna.get_name() << " ---------" << std::endl;

        // Builds loop tree
        knotergy::LoopFactory factory(processed_rna);

        // Computes the energy
        knotergy::ComputeEnergy energy_calculator(factory.get_root_node(), processed_rna, vp, pkp,
                                                  mp, round, verbose);

        // Output results

        if (energy_calculator.getInfiniteEnergyFlag()) {
            std::cout << "\nENERGY: Infinite (" << energy_calculator.getEnergy() << " kcal/mol)"
                      << std::endl;
        } else {
            printf("\nENERGY: %.4f kcal/mol\n", energy_calculator.getEnergy());
        }
    }
    return 0;
}