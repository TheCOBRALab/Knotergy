#include <filesystem>
#include <iostream>
#include <string>

#include "pipeline/input_pipeline.hpp"
#include "preprocessing/RNAEntry.hpp"
#include "preprocessing/RNAProcessedEntry.hpp"

namespace {
void help() {
    std::cout << "Usage: ./Knotergy [options]\n"
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
    knotergy::trim(value);
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
        knotergy::trim(sequence);
    }

    if (!sequence.empty() && !knotergy::validate_sequence(sequence)) {
        std::cout << "Error: Sequence is empty or contains invalid character/s. Allowed: G, C, "
                     "A, U, T";
        return 1;
    }

    if (structure.empty() && input_file.empty()) {
        std::cout << "Structure: ";
        std::cin >> structure;
        knotergy::trim(structure);
    }

    if (!structure.empty() && !knotergy::validate_structure(structure)) {
        std::cout << "Error: Structure is empty or contains invalid character/s. Allowed: '.', "
                     "'(',  ')', '[', ']'";
        return 1;
    }

    if (sequence.length() != structure.length()) {
        std::cout << "Error: Input sequence and structure are not the same length";
        std::cout << "\nSequence length: " << sequence.length()
                  << ", Structure length: " << structure.length() << std::endl;
        return 1;
    }

    knotergy::trim(input_file);
    if (!input_file.empty() && !std::filesystem::exists(input_file)) {
        std::cerr << "Input file not found: " << input_file << std::endl;
        return 1;
    }

    if (!parameter_file.empty() && !std::filesystem::exists(parameter_file)) {
        std::cerr << "Parameter file not found: " << parameter_file << std::endl;
        return 1;
    }

    //------------------------- Pre-processing and reading from files -----------------------------
    std::vector<knotergy::RNAEntry> inputs =
        knotergy::get_all_inputs(input_file, sequence, structure);
    std::vector<knotergy::RNAProcessedEntry> processed_inputs = knotergy::process_inputs(inputs);

    for (const knotergy::RNAProcessedEntry& current : processed_inputs) {
        std::cout << "Name: " << current.get_name() << " Sequence: " << current.get_sequence()
                  << "\nStructure: " << current.get_structure() << std::endl;

        for (size_t n : current.get_pairings()) {
            if (n == knotergy::NULL_INDEX) {
                std::cout << -1 << " ";
            } else {
                std::cout << n << " ";
            }
        }
        std::cout << std::endl;
        knotergy::dostuff(current, parameter_file);
    }
    return 0;
}