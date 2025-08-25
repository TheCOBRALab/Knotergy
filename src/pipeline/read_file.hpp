#pragma once

#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include "shared.hpp"
#include "../preprocessing/RNAEntry.hpp"

namespace knotergy{

/**
 * @brief Represents the parser's state during file processing.
 *
 * Used internally in get_all_file_entries() to track whether the parser
 * is currently reading a name, sequence, or structure.
 */
enum class ParserState { UNINITIALIZED, NAME, SEQUENCE, STRUCTURE };

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
[[nodiscard]] std::vector<RNAEntry> get_all_file_entries(const std::string& file) {
    if (!std::filesystem::exists(file)) {
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
                    THROW_ERROR(
                        "Error: Sequence and/or structure are empty for entry: " + current.name +
                        ". Line number: " + std::to_string(line_number));
                }
                if (current.sequence.size() != current.structure.size()) {
                    THROW_ERROR(
                        "Error: Sequence and structure are not the same length in entry: " +
                        current.name + ". Line number: " + std::to_string(line_number));
                }
                entries.push_back(current);
                current = {};
            }
            current.name = line.substr(1);
            state = ParserState::SEQUENCE;
        } else if (state == ParserState::SEQUENCE) {
            current.sequence = line;
            state = ParserState::STRUCTURE;
        } else if (state == ParserState::STRUCTURE) {
            current.structure = line;
            state = ParserState::NAME;
        } else {
            // Should never reach here
            THROW_ERROR("Error: Unexpected state. Line number: " +
                                     std::to_string(line_number));
        }
    }

    // Saves the last entry
    if (!current.name.empty() && !current.name.empty() && !current.structure.empty()) {
        entries.push_back(current);
    } else if (!current.name.empty() && (current.sequence.empty() || current.structure.empty())) {
        THROW_ERROR("Error: Sequence and/or structure are empty for entry: " +
                                 current.name + ". Line number: " + std::to_string(line_number));
    }

    return entries;
}

}