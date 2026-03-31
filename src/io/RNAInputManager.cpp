#include "RNAInputManager.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../energy/ComputeEnergy.hpp"
#include "../loop_tree/LoopFactory.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "../preprocessing/RNAEntry.hpp"
#include "FileUtils.hpp"
#include "common.hpp"

namespace knotergy {

// Collects all RNA input entries from console and/or file.
std::vector<RNAEntry> RNAInputManager::get_all_inputs(const std::string& input_file,
                                                      const std::string& sequence,
                                                      const std::string& structure) {
    std::vector<RNAEntry> entries;

    // get console input
    if (!sequence.empty()) {
        entries.emplace_back("Console Sequence", sequence, structure);
    }

    // get file inputs
    if (!input_file.empty()) {
        std::vector<RNAEntry> file_entries = RNAInputManager::get_all_FASTA_entries(input_file);
        // move values into entries (avoids deep copies). Keeps console as first entry
        entries.insert(entries.end(), std::make_move_iterator(file_entries.begin()),
                       std::make_move_iterator(file_entries.end()));
    }
    if (entries.empty()) THROW_ERROR("No Input Data Given");
    return entries;
}

std::vector<ProcessedRNAEntry> RNAInputManager::process_inputs(
    const std::vector<RNAEntry>& inputs, const std::vector<modified_base_param>& modified_params) {
    std::vector<ProcessedRNAEntry> processed_inputs;
    processed_inputs.reserve(inputs.size());

    for (const RNAEntry& rna : inputs) {
        processed_inputs.emplace_back(RNAProcessor::process_rna(rna, modified_params));
    }
    return processed_inputs;
}

// Parses RNA entries from a file in FASTA format
[[nodiscard]] std::vector<RNAEntry> RNAInputManager::get_all_FASTA_entries(
    const std::string& file) {
    // Represents the parser's state during file processing.
    // Used to track if it's currently reading name, sequence, or structure.
    enum class ParserState { UNINITIALIZED, NAME, SEQUENCE, STRUCTURE };

    if (!FileUtils::file_exists(file)) {
        THROW_ERROR("Error: Input file not found: " + file);
    }

    std::ifstream in(file);
    if (!in.is_open()) {
        THROW_ERROR("Error: Unable to open file: " + file);
    }

    std::string line;
    RNAEntry current;
    std::vector<RNAEntry> entries;
    ParserState state = ParserState::UNINITIALIZED;
    int line_number = 0;

    // Loop through each line in the file, and is a state machine to parse entries
    // Each entry consists of three lines: name, sequence, structure
    // Each state expects a specific line format (e.g., name line starts with '>')
    while (std::getline(in, line)) {
        ++line_number;
        trim(line);
        if (line.empty()) continue;

        // Name and Uninitialized must start with a header line indicator
        if ((state == ParserState::NAME || state == ParserState::UNINITIALIZED) &&
            (line[0] != '>')) {
            THROW_ERROR("Error: Expected '>' at the beginning of the line: " + line +
                        ". Line number: " + std::to_string(line_number));
        }

        // '>' is a header line indicator, indicating the start of a new entry
        if (line[0] == '>') {
            if (state != ParserState::UNINITIALIZED) {
                if (current.sequence.empty() || current.structure.empty()) {
                    THROW_ERROR("Error: Sequence and/or structure are empty for entry: " +
                                current.name + ". Line number: " + std::to_string(line_number));
                }
                if (current.sequence.size() != current.structure.size()) {
                    THROW_ERROR("Error: Sequence and structure are not the same length in entry: " +
                                current.name + ". Line number: " + std::to_string(line_number));
                }
                entries.push_back(current);
                current = {};
            }
            current.name = line.substr(1);
            state = ParserState::SEQUENCE;
        } else if (state == ParserState::SEQUENCE) {
            current.sequence = std::move(line);
            state = ParserState::STRUCTURE;
        } else if (state == ParserState::STRUCTURE) {
            current.structure = std::move(line);
            state = ParserState::NAME;
        } else {
            // Should never reach here
            THROW_ERROR("Error: Unexpected state. Line number: " + std::to_string(line_number));
        }
    }

    // Saves the last entry
    if (!current.name.empty() && !current.name.empty() && !current.structure.empty()) {
        entries.push_back(current);
    } else if (!current.name.empty() && (current.sequence.empty() || current.structure.empty())) {
        THROW_ERROR("Error: Sequence and/or structure are empty for entry: " + current.name +
                    ". Line number: " + std::to_string(line_number));
    }

    return entries;
}

}  // namespace knotergy