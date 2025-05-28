#include <filesystem>
#include <iostream>
#include <string>

#include "ComputeEnergy.hpp"

namespace {
void help() {
    std::cout << "Usage: ./ComputeEnergy [options]\n"
              << "Options:\n"
              << "  -s, --sequence <string>       RNA sequence\n"
              << "  -r, --structure <string>      Input structure\n"
              << "  -i, --input <file>            Input file\n"
              << "  -o, --output <file>           Output file\n"
              << "  -p, --paramFile <file>        Parameter file\n"
              << "  -h, --help                    Show this help message\n";
}

std::string get_trimmed_arg(int& i, int argc, char** argv) {
    if (i + 1 >= argc) return "";
    std::string value = argv[++i];
    ComputeEnergy::trim(value);
    return value;
}
}  // namespace

int main(int argc, char** argv) {
    std::string sequence = "";
    std::string structure = "";
    std::string input_file = "";
    std::string output_file = "";
    std::string parameter_file = "";

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
        } else if (arg == "-h" || arg == "--help") {
            help();
            return 0;
        } else {
            std::cerr << "Unknown option or missing value: " << arg << std::endl;
            return 1;
        }
    }

    // ------------------------- Validate Inputs -----------------------
    if (sequence.empty() && input_file.empty()) {
        std::cout << "Sequence : ";
        std::cin >> sequence;
        ComputeEnergy::trim(sequence);
    }

    if (!sequence.empty() && !ComputeEnergy::validate_sequence(sequence)) {
        std::cout << "Error: Sequence is empty or contains invalid character/s. Allowed: G, C, "
                     "A, U, T";
        return 1;
    }

    if (structure.empty() && input_file.empty()) {
        std::cout << "Structure: ";
        std::cin >> structure;
        ComputeEnergy::trim(structure);
    }

    if (!structure.empty() && !ComputeEnergy::validate_structure(structure)) {
        std::cout << "Error: Structure is empty or contains invalid character/s. Allowed: '.', "
                     "'(',  ')', '[', ']'";
        return 1;
    }

    if (sequence.length() != structure.length()) {
        std::cout << "Error: Input sequence and structure are not the same length";
        return 1;
    }

    ComputeEnergy::trim(input_file);
    if (!input_file.empty() && !std::filesystem::exists(input_file)) {
        std::cerr << "Input file not found: " << input_file << std::endl;
        return 1;
    }

    if (!parameter_file.empty() && !std::filesystem::exists(parameter_file)) {
        std::cerr << "Parameter file not found: " << parameter_file << std::endl;
        return 1;
    }

    //------------------------- Pre-processing and reading from files -----------------------------
    std::vector<ComputeEnergy::RNAEntry> inputs =
        ComputeEnergy::get_all_inputs(input_file, sequence, structure);

    for (ComputeEnergy::RNAEntry& current : inputs) {
        // std::cout << "Name: " << current.name << " Sequence: " << current.sequence << "
        // Structure: " << current.structure << std::endl;
    }
    return 0;
}