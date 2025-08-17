#include "input_pipeline.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../energy/ComputeEnergy.hpp"
#include "../loop_tree/LoopFactory.hpp"
#include "../preprocessing/RNAEntry.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "shared.hpp"



namespace knotergy {

// Trims leading and trailing whitespace from a string
void trim(std::string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
    s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
}

// Gets all inputs from a file in FASTA format
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
                if (current.sequence.empty() || current.structure.empty()) {
                    throw std::runtime_error(
                        "Error: Sequence and/or structure are empty for entry: " + current.name +
                        ". Line number: " + std::to_string(line_number));
                }
                if (current.sequence.size() != current.structure.size()) {
                    throw std::runtime_error(
                        "Error: Sequence and structure are not the same length in entry: " +
                        current.name + ". Line number: " + std::to_string(line_number));
                }
                entries.push_back(current);
                current = {};
            }
            current.name = line.substr(1);
            state = ParserState::SEQUENCE;
        } else if (state == ParserState::SEQUENCE) {
            if (!validate_sequence(line)) {
                throw std::runtime_error("Error: Sequence is invalid for entry: " + current.name +
                                         ". Line number: " + std::to_string(line_number));
            }
            current.sequence = line;
            state = ParserState::STRUCTURE;
        } else if (state == ParserState::STRUCTURE) {
            // if (!validate_structure(line)) {
            //     throw std::runtime_error("Error: Structure is invalid for entry: " + current.name
            //     +
            //                              ". Line number: " + std::to_string(line_number));
            // }
            current.structure = line;
            state = ParserState::NAME;
        } else {
            // Should never reach here
            throw std::runtime_error("Error: Unexpected state. Line number: " +
                                     std::to_string(line_number));
        }
    }

    // Saves the last entry
    if (!current.name.empty() && !current.name.empty() && !current.structure.empty()) {
        entries.push_back(current);
    } else if (!current.name.empty() && (current.sequence.empty() || current.structure.empty())) {
        throw std::runtime_error("Error: Sequence and/or structure are empty for entry: " +
                                 current.name + ". Line number: " + std::to_string(line_number));
    }

    return entries;
}

// Collects all RNA input entries from console input and/or file.
std::vector<RNAEntry> get_all_inputs(const std::string& input_file, const std::string& sequence,
                                     const std::string& structure) {
    std::vector<RNAEntry> entries;

    // get console input
    if (!sequence.empty()) {
        if (sequence.size() != structure.size()) {
            THROW_ERROR("Input sequence and structure are not the same length.\nSequence length: " +
                        std::to_string(sequence.size()) +
                        "\nStructure length: " + std::to_string(structure.size()));
        }
        entries.emplace_back("Console Sequence", sequence, structure);
    }

    // get file inputs
    if (!input_file.empty()) {
        std::vector<RNAEntry> file_entries = get_all_file_entries(input_file);
        // move values into entries (avoids deep copies). Keeps console as first entry
        entries.insert(entries.end(), std::make_move_iterator(file_entries.begin()),
                       std::make_move_iterator(file_entries.end()));
    }
    if (entries.empty()) THROW_ERROR("No Input Data Given");
    return entries;
}

std::vector<ProcessedRNAEntry> process_inputs(const std::vector<RNAEntry>& inputs) {
    std::vector<ProcessedRNAEntry> processed_inputs;
    processed_inputs.reserve(inputs.size());
    for (RNAEntry rna : inputs) {
        processed_inputs.emplace_back(RNAProcessor::process_rna(std::move(rna)));
    }
    return processed_inputs;
}

// ensures sequence only has valid characters
bool validate_sequence(const std::string& sequence) {
    for (char c : sequence) {
        if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
            return false;
        }
    }
    return !sequence.empty();
}

// ensures structure is balanced
[[nodiscard]] bool validate_structure(const std::string& structure, bool throw_error) {
    if (structure.empty()) {
        if (throw_error) {
            THROW_ERROR("Invalid RNA structure: Structure is empty");
        }
        return false;
    }

    std::unordered_map<char, char> open_to_close = {{'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}};

    std::unordered_map<char, char> close_to_open = {{')', '('}, {']', '['}, {'}', '{'}, {'>', '['}};

    std::unordered_map<char, int> open_count = {{'(', 0}, {'[', 0}, {'{', 0}, {'<', 0}};

    for (size_t i = 0; i < structure.size(); i++) {
        char c = structure[i];

        if (c == '.') continue;

        if (open_to_close.count(c)) {
            ++open_count[c];
            continue;
        }

        if (close_to_open.count(c)) {
            char open = close_to_open[c];
            if (open_count[open] <= 0) {
                if (throw_error) {
                    THROW_ERROR("Invalid RNA structure: Bracket: '" + std::string(1, c) +
                                "' at index: " + std::to_string(i) + " was never opened");
                }
                return false;
            }
            --open_count[c];
        }
    }

    for (auto& [open, count] : open_count) {
        if (count <= 0) {
            if (throw_error) {
                THROW_ERROR("Invalid RNA structure: opening bracket '" + std::string(1, open) +
                            "' was not closed");
            }
            return false;
        }
    }
    return true;
}

void dostuff(const ProcessedRNAEntry& processed_rna, std::string parameter_file, bool round) {
    ViennaParams::load_energy_parameters(parameter_file, processed_rna.get_sequence());
    printf("Seq: %s \n", processed_rna.get_sequence().c_str());
    printf("Struct: %s \n", processed_rna.get_structure().c_str());
    printf("Size: %ld \n", processed_rna.size());
    for (size_t i = 1; i < processed_rna.size(); ++i) {
        // printf("%d ", entry.get_pairings()[i]);
    }
    printf("\n-------------------------------\n Making the Loop Tree\n");
    LoopFactory factory(processed_rna);
    factory.print_tree(true);
    ComputeEnergy energy_calculator(factory.get_root_node(), processed_rna.get_sequence(),
                                    processed_rna, round);
    std::cout << "ENERGY: " << energy_calculator.getEnergy() << std::endl;
}
}  // namespace knotergy