#ifndef COMPUTEENERGY_COMPUTEENERGY_H
#define COMPUTEENERGY_COMPUTEENERGY_H

#include <string>
#include <vector>

namespace ComputeEnergy {
struct RNAEntry {
    std::string name;
    std::string sequence;
    std::string structure;

    // Constructor that takes all three fields
    RNAEntry(std::string n, std::string s, std::string st)
        : name(std::move(n)), sequence(std::move(s)), structure(std::move(st)) {}

    // Default constructor (needed for vector resizing or default initialization)
    RNAEntry() = default;
};

void trim(std::string& s);
[[nodiscard]] bool validate_sequence(const std::string& sequence);
[[nodiscard]] bool validate_structure(const std::string& structure);
std::vector<RNAEntry> get_all_file_entries(const std::string& file);
std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq,
                                     const std::string& restricted);
std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq,
                                     const std::string& restricted);
}  // namespace ComputeEnergy

#endif  // COMPUTEENERGY_COMPUTEENERGY_H