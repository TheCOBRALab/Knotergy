#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "../preprocessing/RNAEntry.hpp"
#include "../preprocessing/RNAProcessor.hpp"
#include "ViennaParams.hpp"
#include "common.hpp"

namespace knotergy {

/**
 * @brief Manages input of RNA sequences and structures from various sources.
 *
 * This class provides static methods to read RNA data from files (FASTA format)
 * or command-line inputs, and to process them into ProcessedRNAEntry objects.
 */
class RNAInputManager {
   public:
    /**
     * @brief Collects all RNA input entries from console input and/or file.
     *
     * This function combines inputs from both user-supplied sequences and
     * file-based entries. If both are provided, all are included.
     *
     * @param fileI Path to input file (optional, can be empty).
     * @param seq A raw RNA sequence string (optional, can be empty).
     * @param restricted Optional structure string associated with the raw sequence.
     * @return A vector of all collected RNAEntry objects.
     *
     * @throws std::runtime_error if both the file and console inputs are empty.
     */
    [[nodiscard]] static std::vector<RNAEntry> get_all_inputs(const std::string& fileI,
                                                              const std::string& seq,
                                                              const std::string& restricted);

    /**
     * @brief Process RNA entries into ProcessedRNAEntry objects.
     *
     * Converts RNAEntry objects into ProcessedRNAEntry objects with precomputed
     * structural annotations (pairings, closed regions, unpaired counts).
     *
     * @param inputs Vector of RNAEntry objects to process.
     * @param modified_params Vector of modified base parameters (default: empty).
     * @return Vector of ProcessedRNAEntry objects.
     */
    [[nodiscard]] static std::vector<ProcessedRNAEntry> process_inputs(
        const std::vector<RNAEntry>& inputs,
        const std::vector<modified_base_param>& modified_params = {});

   private:
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
    [[nodiscard]] static std::vector<RNAEntry> get_all_FASTA_entries(const std::string& file);
};
}  // namespace knotergy
