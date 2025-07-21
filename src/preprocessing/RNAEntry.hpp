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
    RNAEntry(std::string name, std::string sequence, std::string structure)
        : name{std::move(name)}, sequence{std::move(sequence)}, structure{std::move(structure)} {}

    // Constructor without name (defaults to "N/A")
    RNAEntry(std::string sequence, std::string structure)
        : RNAEntry("N/A", std::move(sequence), std::move(structure)) {}

    size_t size() const {
        assert(sequence.size() == structure.size() &&
               "Sequence and structure must be the same length");
        return structure.size();
    }
};

}  // namespace knotergy
