#pragma once

#include <cassert>
#include <string>

namespace knotergy {

/**
 * @brief Represents an RNA entry with name, sequence, and structure.
 *
 * This struct is designed to store information about an RNA molecule, 
 * including its name, nucleotide sequence, and secondary structure. 
 * It also provides utility functions for checking the size of the RNA entry.
 */

struct RNAEntry {
    std::string name;
    std::string sequence;
    std::string structure;

    RNAEntry() = default;

    RNAEntry(std::string rna_name, std::string rna_sequence, std::string rna_structure)
        : name{rna_name},
          sequence{rna_sequence},
          structure{rna_structure} {}

    RNAEntry(std::string rna_sequence, std::string rna_structure)
        : RNAEntry("N/A", rna_sequence, rna_structure) {}

    size_t size() const {
        return structure.size();
    }
};

}  // namespace knotergy
