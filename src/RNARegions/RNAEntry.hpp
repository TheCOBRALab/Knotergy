#pragma once

#include <iostream>
#include <stack>
#include <string>
#include <vector>

namespace ComputeEnergy {
class RNAEntry {
   public:
    RNAEntry(std::string name_, std::string sequence_, std::string structure_) {
        name = std::move(name_);
        sequence = std::move(sequence_);
        set_structure(std::move(structure_));  // makes call to update pairings
    }

    // Default constructor (needed for vector resizing or default initialization)
    RNAEntry() = default;

    // Getters
    [[nodiscard]] const std::string& get_name() const { return name; }
    [[nodiscard]] const std::string& get_sequence() const { return sequence; }
    [[nodiscard]] const std::string& get_structure() const { return structure; }
    [[nodiscard]] size_t get_first_pair_index() const { return first_pair_index; }
    [[nodiscard]] const std::vector<int>& get_pairings() const {
        if (sequence.size() != structure.size()) {
            std::cerr << "Warning: Sequence and Structure are different sizes.\n"
                      << "Sequence: " << sequence << "\n"
                      << "Structure: " << structure << std::endl;
        }
        return pairings;
    }

    // Setters
    void set_name(std::string name_) { name = name_; }
    void set_sequence(std::string sequence_) { sequence = sequence_; }
    void set_structure(std::string structure_) {
        structure = structure_;
        update_pairings();
    }

   private:
    std::string name;
    std::string sequence;
    std::string structure;
    std::vector<int> pairings;  // [4, -1, -1, 1] Number represents the index that base is paired to
    size_t first_pair_index = 0;

    /**
     * @brief Computes base pairings from the RNA secondary structure string.
     *
     * This function parses the `structure` string using dot-bracket notation and updates
     * the `pairings` vector, where each index `i` maps to the base `j` that `i` is paired with.
     * Unpaired positions are marked with `-1`.
     *
     * Supported brackets:
     * - `()` for regular base pairs
     * - `[]` for pseudoknots
     *
     * If the `structure` contains mismatched brackets, unrecognized characters, or is
     * inconsistent with the length of `sequence`, this function throws a `std::runtime_error`.
     *
     * Warnings:
     * - If `sequence` and `structure` differ in length, a warning is printed to `std::cerr`.
     *
     * Exceptions:
     * - `std::runtime_error` if:
     *   - Closing brackets have no matching opener
     *   - Invalid characters are encountered
     *   - Opening brackets remain unmatched at the end
     */
    void update_pairings() {
        pairings.assign(structure.size(), -1);  // pre-allocate pairings with -1
        std::stack<size_t> brackets;
        std::stack<size_t> pseudoknots;
        size_t j;

        for (size_t i = 0; i < structure.size(); i++) {
            switch (structure[i]) {
                case '.':
                    break;
                case '(':
                    if (!first_pair_index) {
                        first_pair_index = i;
                    }
                    brackets.push(i);
                    break;
                case '[':
                    if (!first_pair_index) {
                        first_pair_index = i;
                    }
                    pseudoknots.push(i);
                    break;
                case ')':
                    if (brackets.empty()) {
                        throw std::runtime_error("Structure in RNAEntry is invalid. \nSequence: " +
                                                 sequence + "\nStructure: " + structure);
                    }
                    j = brackets.top();
                    brackets.pop();
                    pairings[i] = static_cast<int>(j);
                    pairings[j] = static_cast<int>(i);
                    break;
                case ']':
                    if (pseudoknots.empty()) {
                        throw std::runtime_error("Structure in RNAEntry is invalid. \nSequence: " +
                                                 sequence + "\nStructure: " + structure);
                    }
                    j = pseudoknots.top();
                    pseudoknots.pop();
                    pairings[i] = static_cast<int>(j);
                    pairings[j] = static_cast<int>(i);
                    break;
                default:
                    throw std::runtime_error(
                        "Character in RNAEntry's structure is invalid. \nInvalid Character: " +
                        std::string(1, structure[i]) + "\nSequence: " + sequence +
                        "\nStructure: " + structure);
            }
        }
        if (!brackets.empty() || !pseudoknots.empty()) {
            throw std::runtime_error("Unmatched opening brackets in RNA structure.\nSequence: " +
                                     sequence + "\nStructure: " + structure);
        }
    }
};
}  // namespace ComputeEnergy