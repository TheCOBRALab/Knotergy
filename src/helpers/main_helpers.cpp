#include "main_helpers.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../energy/ComputeEnergy.hpp"
#include "../loops/LoopFactory.hpp"
#include "../rna_regions/RNAEntry.hpp"
#include "common.hpp"

extern "C" {
// used for load_energy_parameters
#include <ViennaRNA/params/io.h>
}

namespace knotergy {

/**
 * @brief Trims leading and trailing whitespace from a string.
 *
 * This function modifies the input string in-place to remove any leading
 * and trailing whitespace characters, including spaces, tabs, newlines,
 * carriage returns, form feeds, and vertical tabs.
 *
 * @param s The string to be trimmed.
 */
void trim(std::string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
    s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
}

/**
 * @brief Parses RNA entries from a file in FASTA format
 *
 * Each entry is expected to consist of three lines in the following order:
 * - Name line starting with '>'
 * - RNA sequence line (ACGTU only)
 * - Structure line (dot-bracket notation)
 *
 * If any line is malformed or a required field is missing, an exception is thrown.
 *
 * @param file Path to the input file containing RNA entries.
 * @return A vector of RNAEntry objects parsed from the file.
 *
 * @throws std::runtime_error if the file does not exist, is unreadable,
 *         or contains malformed data.
 */
std::vector<RNAEntry> get_all_file_entries(const std::string& file) {
    if (!std::filesystem::exists(file)) {
        throw std::runtime_error("Error: Input file not found: " + file);
    }

    std::ifstream in(file);
    if (!in.is_open()) {
        throw std::runtime_error("Error: Unable to open file: " + file);
    }

    std::string line;
    RNAEntry current;
    std::vector<RNAEntry> entries;
    ParserState state = ParserState::UNINITIALIZED;
    int line_number = 0;

    while (std::getline(in, line)) {
        ++line_number;
        trim(line);
        if (line.empty()) continue;

        // Name and Uninitialized must start with a header line indicator
        if ((state == ParserState::NAME || state == ParserState::UNINITIALIZED) &&
            (line[0] != '>')) {
            throw std::runtime_error("Error: Expected '>' at the beginning of the line: " + line +
                                     ". Line number: " + std::to_string(line_number));
        }

        // '>' is a header line indicator, indicating the start of a new entry
        if (line[0] == '>') {
            if (state != ParserState::UNINITIALIZED) {
                if (current.get_sequence().empty() || current.get_structure().empty()) {
                    throw std::runtime_error(
                        "Error: Sequence and/or structure are empty for entry: " +
                        current.get_name() + ". Line number: " + std::to_string(line_number));
                }
                if (current.get_sequence().length() != current.get_structure().length()) {
                    throw std::runtime_error(
                        "Error: Sequence and structure are not the same length in entry: " +
                        current.get_name() + ". Line number: " + std::to_string(line_number));
                }
                entries.push_back(current);
                current = {};
            }
            current.set_name(line.substr(1));
            state = ParserState::SEQUENCE;
        } else if (state == ParserState::SEQUENCE) {
            if (!validate_sequence(line)) {
                throw std::runtime_error(
                    "Error: Sequence is invalid for entry: " + current.get_name() +
                    ". Line number: " + std::to_string(line_number));
            }
            current.set_sequence(line);
            state = ParserState::STRUCTURE;
        } else if (state == ParserState::STRUCTURE) {
            if (!validate_structure(line)) {
                throw std::runtime_error(
                    "Error: Structure is invalid for entry: " + current.get_name() +
                    ". Line number: " + std::to_string(line_number));
            }
            current.set_structure(line);
            state = ParserState::NAME;
        } else {
            // Should never reach here
            throw std::runtime_error("Error: Unexpected state. Line number: " +
                                     std::to_string(line_number));
        }
    }

    // Saves the last entry
    if (!current.get_name().empty() && !current.get_sequence().empty() &&
        !current.get_structure().empty()) {
        entries.push_back(current);
    } else if (!current.get_name().empty() &&
               (current.get_sequence().empty() || current.get_structure().empty())) {
        throw std::runtime_error(
            "Error: Sequence and/or structure are empty for entry: " + current.get_name() +
            ". Line number: " + std::to_string(line_number));
    }

    return entries;
}

/**
 * @brief Collects all RNA input entries from console input and/or file.
 *
 * This function combines inputs from both user-supplied sequences and
 * file-based entries. If both are provided, all are included.
 *
 * @param input_file Path to input file (optional, can be empty).
 * @param sequence A raw RNA sequence string (optional, can be empty).
 * @param structure Optional structure string associated with the raw sequence.
 * @return A vector of all collected RNAEntry objects.
 *
 * @throws std::runtime_error if both the file and console inputs are empty
 */
std::vector<RNAEntry> get_all_inputs(const std::string& input_file, const std::string& sequence,
                                     const std::string& structure) {
    std::vector<RNAEntry> entries;
    if (!sequence.empty()) {
        entries.emplace_back("Console Sequence", sequence, structure);
    }
    if (!input_file.empty()) {
        std::vector<RNAEntry> file_entries = get_all_file_entries(input_file);
        entries.insert(entries.end(), std::make_move_iterator(file_entries.begin()),
                       std::make_move_iterator(file_entries.end()));
    }
    if (entries.empty()) throw std::runtime_error("No Input Data Given");
    return entries;
}

/**
 * @brief Validates that an RNA sequence contains only valid characters.
 *
 * Accepted characters are A, C, G, U, and T.
 *
 * @param sequence The RNA sequence to validate.
 * @return true if the sequence is valid and non-empty; false otherwise.
 */
[[nodiscard]] bool validate_sequence(const std::string& sequence) {
    for (char c : sequence) {
        if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
            return false;
        }
    }
    return !sequence.empty();
}

/**
 * @brief Validates the structure string for RNA using dot-bracket notation.
 *
 * The function ensures that parentheses and square brackets are correctly matched and balanced.
 * Accepted characters: '.', '(', ')', '[', ']'.
 *
 * @param structure The structure string to validate.
 * @return true if the structure is valid and non-empty; false otherwise.
 */
[[nodiscard]] bool validate_structure(const std::string& structure) {
    int paren = 0, square = 0;
    for (char c : structure) {
        switch (c) {
            case '.':
                break;
            case '(':
                ++paren;
                break;
            case ')':
                if (--paren < 0) return false;
                break;
            case '[':
                ++square;
                break;
            case ']':
                if (--square < 0) return false;
                break;
            default:
                return false;
        }
    }
    return !structure.empty() && paren == 0 && square == 0;
}

void load_energy_parameters(const std::string& paramFile, const std::string& seq) {
    if (!paramFile.empty()) {
        if (std::filesystem::exists(paramFile)) {
            int loaded = vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            if (!loaded) {
                throw std::runtime_error("Failed to load parameter file: " + paramFile);
            }
            std::cout << "Successfully loaded parameter file: " << paramFile << std::endl;
            return;
        } else {
            std::cerr << "Warning: Parameter file \"" << paramFile << "\" not found." << std::endl;
        }
    } else {
        std::cerr << "Warning: No parameter file provided." << std::endl;
    }

    // Default fallback based on sequence
    if (seq.find('T') != std::string::npos) {
        std::cerr << "Defaulting to DNA parameters (Mathews 2004)." << std::endl;
        vrna_params_load_DNA_Mathews2004();
    } else {
        std::cerr << "Defaulting to RNA parameters (Langdon 2018)." << std::endl;
        vrna_params_load_RNA_Langdon2018();
    }
}

void dostuff(RNAEntry entry, std::string parameter_file) {
    load_energy_parameters(parameter_file, entry.get_sequence());
    printf("Seq: %s \n", entry.get_sequence().c_str());
    printf("Size: %ld \n", entry.get_sequence().size());
    for (size_t i = 1; i < entry.get_sequence().size(); ++i) {
        // printf("%d ", entry.get_pairings()[i]);
    }
    printf("\n-------------------------------\n Making the Loop Tree\n");
    LoopFactory factory(entry);
    // factory.print_tree(true);
    ComputeEnergy energy_calculator(factory.get_root_node(), entry.get_sequence());
    std::cout << "ENERGY: " << energy_calculator.getEnergy() << std::endl;
}

}  // namespace knotergy