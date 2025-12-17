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

[[nodiscard]] std::vector<ProcessedRNAEntry> process_inputs(const std::vector<RNAEntry>& inputs,
                                 const std::vector<modified_base_params>& modified_params = {});

}  // namespace knotergy
