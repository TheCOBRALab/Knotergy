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
    std::string name;       ///< Name or identifier of the RNA entry.
    std::string sequence;   ///< RNA nucleotide sequence.
    std::string structure;  ///< Dot-bracket structure notation.

    RNAEntry() = default;

    /**
     * @brief Construct an RNA entry with name, sequence, and structure.
     *
     * @param rna_name Name or identifier for the RNA.
     * @param rna_sequence Nucleotide sequence.
     * @param rna_structure Dot-bracket structure string.
     */
    RNAEntry(std::string rna_name, std::string rna_sequence, std::string rna_structure)
        : name{rna_name}, sequence{rna_sequence}, structure{rna_structure} {}

    /**
     * @brief Construct an RNA entry with sequence and structure (name defaults to "N/A").
     *
     * @param rna_sequence Nucleotide sequence.
     * @param rna_structure Dot-bracket structure string.
     */
    RNAEntry(std::string rna_sequence, std::string rna_structure)
        : RNAEntry("N/A", rna_sequence, rna_structure) {}

    /**
     * @brief Get the length of the RNA structure.
     *
     * @return Size of the structure string.
     */
    size_t size() const { return structure.size(); }
};

}  // namespace knotergy
