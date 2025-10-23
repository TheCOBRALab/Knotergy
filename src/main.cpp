#include <iostream>
#include <string>
#include <unordered_set>

#include "pipeline/input_pipeline.hpp"
#include "preprocessing/RNAEntry.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
#include "loop_tree/LoopFactory.hpp"
#include "energy/ComputeEnergy.hpp"

#define KNOTERGY_VERSION "0.1.2"

namespace {
void help() {
    std::cout << "Usage: ./Knotergy [options]\n"
              << "Options:\n"
              << "  -h, --help                    Show this help message\n"
              << "  -V, --version                 Print version and exit\n"
              << "  -s, --sequence <string>       RNA sequence\n"
              << "  -r, --structure <string>      Input structure\n"
              << "  -i, --input <file>            Input file\n"
              << "  -o, --output <file>           Output file\n"
              << "  -p, --paramFile <file>        Parameter file\n"
              << "  -e  --round                   Rounds all decimal places in pseudoknot calculations";
            //   << "  -m, --modifications           Chars used for modified bases (default=`7I6P9D')\n"
            //   << "  -f, --mod-file                Modified base parameter file\n"
            //   << "  -d, --dangle                  Specify the dangle model to be used (base is 2)\n";

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
    std::string mod_param_file = "";
    bool round = false;
    int  dangle = 2;

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
        } else if (arg == "-m" || arg == "--modifications") {
            std::cerr << "Modified bases are not currently supported. This flag will be ignored." << std::endl;
            modifications = get_trimmed_arg(i, argc, argv);
        } else if (arg == "-f" || arg == "--mod-file") {
            std::cerr << "Modified bases are not currently supported. This flag will be ignored." << std::endl;
            mod_param_file = get_trimmed_arg(i, argc, argv);
        } else if (arg == "-d" || arg == "--dangle") {
            std::cerr << "Dangles are not currently supported. This flag will be ignored." << std::endl;
            dangle = get_numerical_arg(i, argc, argv, dangle);
        } else if (arg == "-h" || arg == "--help") {
            help();
            return 0;
        } else if (arg == "-v" || arg == "--version") {
            std::cout << "Knotergy " << KNOTERGY_VERSION << std::endl;
            return 0;
        } else {
            std::cerr << "Unknown option or missing value: " << arg << std::endl;
            return 1;
        }
    }

    std::unordered_set<char> valid_seq_chars(modifications.begin(), modifications.end());
    valid_seq_chars.insert({'A', 'U', 'C', 'G', 'T'});

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

    if (sequence.length() != structure.length()) {
        std::cout << "Error: Input sequence and structure are not the same length";
        std::cout << "\nSequence length: " << sequence.length()
                  << ", Structure length: " << structure.length() << std::endl;
        return 1;
    }

    knotergy::trim(input_file);
    if (!input_file.empty() && !knotergy::file_exists(input_file)) {
        std::cerr << "Input file not found: " << input_file << std::endl;
        return 1;
    }

    knotergy::trim(parameter_file);
    if (!parameter_file.empty() && !knotergy::file_exists(parameter_file)) {
        std::cerr << "Parameter file not found: " << parameter_file << std::endl;
        return 1;
    }

    knotergy::trim(mod_param_file);
    if (!mod_param_file.empty() && !knotergy::file_exists(mod_param_file)) {
        std::cerr << "Modified bases parameter file not found: " << mod_param_file << std::endl;
        return 1;
    }


    //------------------------- Pre-processing and reading from files -----------------------------
    std::vector<knotergy::RNAEntry> inputs = knotergy::get_all_inputs(input_file, sequence, structure);
    std::vector<knotergy::ProcessedRNAEntry> processed_inputs = knotergy::process_inputs(inputs);
    knotergy::ViennaParams::load_energy_parameters(parameter_file);

    //------------------------- Main Processing Loop -----------------------------
    for (const knotergy::ProcessedRNAEntry& processed_rna : processed_inputs) {
        knotergy::validate_sequence(sequence, valid_seq_chars);

        knotergy::LoopFactory factory(processed_rna);
        knotergy::ComputeEnergy energy_calculator(factory.get_root_node(), processed_rna.get_sequence(), processed_rna, round);
        // std::cout << "\nName: " << processed_rna.get_name() << "\nSequence: " << processed_rna.get_sequence()
        //           << "\nStructure: " << processed_rna.get_structure() << std::endl;
        printf("\nENERGY: %.4f kcal/mol\n", energy_calculator.getEnergy());
    }
    return 0;
}