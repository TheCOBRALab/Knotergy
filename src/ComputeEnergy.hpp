#pragma once

#include <string>
#include <vector>

#include "RNARegions/RNAEntry.hpp"

namespace ComputeEnergy {

struct Region {
    int begin = -1;
    int end = -1;
};

/**
 * @brief Represents the parser's state during file processing.
 *
 * Used internally in get_all_file_entries() to track whether the parser
 * is currently reading a name, sequence, or structure.
 */
enum class ParserState { UNINITIALIZED = -1, NAME = 0, SEQUENCE = 1, STRUCTURE = 2 };

void trim(std::string& s);
[[nodiscard]] bool validate_sequence(const std::string& sequence);
[[nodiscard]] bool validate_structure(const std::string& structure);
std::vector<RNAEntry> get_all_file_entries(const std::string& file);
std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq,
                                     const std::string& restricted);
void dostuff(RNAEntry entry);
}  // namespace ComputeEnergy
