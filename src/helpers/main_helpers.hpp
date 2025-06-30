#pragma once

#include <string>
#include <vector>

#include "../rna_regions/RNAEntry.hpp"

namespace knotergy {

/**
 * @brief Represents the parser's state during file processing.
 *
 * Used internally in get_all_file_entries() to track whether the parser
 * is currently reading a name, sequence, or structure.
 */
enum class ParserState { UNINITIALIZED, NAME, SEQUENCE, STRUCTURE };

void trim(std::string& s);
[[nodiscard]] bool validate_sequence(const std::string& sequence);
[[nodiscard]] bool validate_structure(const std::string& structure);
std::vector<RNAEntry> get_all_file_entries(const std::string& file);
std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq,
                                     const std::string& restricted);
void dostuff(RNAEntry entry, std::string parameter_file);
}  // namespace knotergy
