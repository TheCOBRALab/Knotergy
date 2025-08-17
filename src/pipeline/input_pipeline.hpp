#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "../preprocessing/RNAEntry.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "../preprocessing/RNAProcessor.hpp"
#include "load_params.hpp"
#include "shared.hpp"

namespace knotergy {

/**
 * @brief Represents the parser's state during file processing.
 *
 * Used internally in get_all_file_entries() to track whether the parser
 * is currently reading a name, sequence, or structure.
 */
enum class ParserState { UNINITIALIZED, NAME, SEQUENCE, STRUCTURE };

/**
 * @brief Trims leading and trailing whitespace from a string.
 *
 * This function modifies the input string in-place to remove any leading
 * and trailing whitespace characters, including spaces, tabs, newlines,
 * carriage returns, form feeds, and vertical tabs.
 *
 * @param s The string to be trimmed.
 */
void trim(std::string& s);

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
[[nodiscard]] std::vector<RNAEntry> get_all_file_entries(const std::string& file);

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
[[nodiscard]] std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq,
                                                   const std::string& restricted);

[[nodiscard]] std::vector<ProcessedRNAEntry> process_inputs(const std::vector<RNAEntry>& inputs);

/**
 * @brief Validates that an RNA sequence contains only valid characters.
 *
 * Accepted characters are A, C, G, U, and T.
 *
 * @param sequence The RNA sequence to validate.
 * @return true if the sequence is valid and non-empty; false otherwise.
 */
[[nodiscard]] bool validate_sequence(const std::string& sequence);

void dostuff(const ProcessedRNAEntry& entry, std::string parameter_file, bool round = false);
}  // namespace knotergy
