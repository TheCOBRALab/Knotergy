#pragma once

#include <string>
#include <vector>

#include "../preprocessing/RNAEntry.hpp"
#include "../preprocessing/RNAProcessedEntry.hpp"
#include "shared.hpp"

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
std::vector<RNAProcessedEntry> process_inputs(const std::vector<RNAEntry>& inputs);

void load_energy_parameters(const std::string& paramFile, const std::string& seq);
void load_energy_parameters(const std::string& paramFile);
void load_energy_parameters();

void dostuff(const RNAProcessedEntry& entry, std::string parameter_file);
}  // namespace knotergy
