#pragma once

#include <iostream>
#include <map>
#include <stack>
#include <string>
#include <vector>

#include "../helpers/common.hpp"

namespace compute_energy {
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
    [[nodiscard]] const std::vector<size_t>& get_pairings() const {
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
        generate_unpaired_bases_count_map();
    }

    // [from, to)
    size_t get_unpaired_count(size_t from, size_t to) {
        if (from >= to) return 0;
        // todo: validate input
        return unpaired_count_map[to] - unpaired_count_map[from];
    }

   private:
    std::string name;
    std::string sequence;
    std::string structure;
    std::vector<size_t>
        pairings;  // [4, -1, -1, 1] Number represents the index that base is paired to
    std::map<size_t, size_t> unpaired_count_map;

    /**
     * @brief Computes base pairings from the RNA secondary structure string.
     *
     * This function loops through a structure in dot-bracket notation
     * It populates the vector `pairings` with the indicies of each pair
     * for example
     *  .(.[.).]..
     *  -1, 5, -1, 7, -1, 1, -1, 3, -1, -1
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
        pairings.assign(structure.size(), NULL_INDEX);  // pre-allocate pairings with -1
        std::stack<size_t> brackets;
        std::stack<size_t> pseudoknots;
        size_t j;

        for (size_t i = 0; i < structure.size(); i++) {
            switch (structure[i]) {
                case '.':
                    break;
                case '(':
                    brackets.push(i);
                    break;
                case '[':
                    pseudoknots.push(i);
                    break;
                case ')':
                    if (brackets.empty()) {
                        throw std::runtime_error("Structure in RNAEntry is invalid. \nSequence: " +
                                                 sequence + "\nStructure: " + structure);
                    }
                    j = brackets.top();
                    brackets.pop();
                    pairings[i] = j;
                    pairings[j] = i;
                    break;
                case ']':
                    if (pseudoknots.empty()) {
                        throw std::runtime_error("Structure in RNAEntry is invalid. \nSequence: " +
                                                 sequence + "\nStructure: " + structure);
                    }
                    j = pseudoknots.top();
                    pseudoknots.pop();
                    pairings[i] = j;
                    pairings[j] = i;
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

    /**
     * @brief Creates a map indicating the number of unpaired bases up till that index
     */
    void generate_unpaired_bases_count_map() {
        size_t count = 0;
        for (size_t i = 0; i <= structure.size(); i++) {
            unpaired_count_map[i] = count;
            if (pairings[i] == NULL_INDEX) {
                ++count;
            }
        }
    };
};
}  // namespace compute_energy