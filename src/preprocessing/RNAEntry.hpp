#pragma once

#include <cassert>
#include <string>

namespace knotergy {

struct RNAEntry {
    std::string name;
    std::string sequence;
    std::string structure;

    // Default constructor (required for containers and default initialization)
    RNAEntry() = default;

    // Constructor with all fields
    RNAEntry(std::string rna_name, std::string rna_sequence, std::string rna_structure)
        : name{std::move(rna_name)}, sequence{std::move(rna_sequence)}, structure{std::move(rna_structure)} {}

    // Constructor without name (defaults to "N/A")
    RNAEntry(std::string rna_sequence, std::string rna_structure)
        : RNAEntry("N/A", std::move(rna_sequence), std::move(rna_structure)) {}

    size_t size() const {
        assert(sequence.size() == structure.size() &&
               "Sequence and structure must be the same length");
        return structure.size();
    }
};

}  // namespace knotergy
