#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../preprocessing/RNAEntry.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "../preprocessing/RNAProcessor.hpp"
#include "load_params.hpp"
#include "shared.hpp"

namespace knotergy {

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
 * Accepted characters are A, C, G, U, T & all modified bases
 *
 * @param sequence The RNA sequence to validate.
 * @param valid_seq_chars Valid chars including modified bases
 * @throws std::runtime_error if the sequence contains invalid characters.
 * @return None
 */
void validate_sequence(const std::string& sequence, const std::unordered_set<char>& valid_seq_chars);

}  // namespace knotergy
