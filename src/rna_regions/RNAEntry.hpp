#pragma once

#include <iostream>
#include <stack>
#include <string>
#include <vector>

#include "../helpers/common.hpp"

namespace compute_energy {
class RNAEntry {
   public:
    RNAEntry(std::string name, std::string sequence, std::string structure) {
        name_ = std::move(name);
        sequence_ = std::move(sequence);
        set_structure(std::move(structure));  // makes call to update pairings
    }

    // Default constructor (needed for vector resizing or default initialization)
    RNAEntry() = default;

    // Getters
    [[nodiscard]] const std::string& get_name() const { return name_; }
    [[nodiscard]] const std::string& get_sequence() const { return sequence_; }
    [[nodiscard]] const std::string& get_structure() const { return structure_; }
    [[nodiscard]] const std::vector<bool>& get_pseudoknot_flags() const {
        return is_pseudoknot_pair;
    }
    [[nodiscard]] const std::vector<Pair>& get_pair_struct() const {
        return pairs_struct;
    }
    [[nodiscard]] const std::vector<size_t>& get_pairings() const {
        if (sequence_.size() != structure_.size()) {
            std::cerr << "Warning: Sequence and Structure are different sizes.\n"
                      << "Sequence: " << sequence_ << "\n"
                      << "Structure: " << structure_ << std::endl;
        }
        return pairings;
    }

    // Setters
    void set_name(std::string name) { name_ = name; }
    void set_sequence(std::string sequence) { sequence_ = sequence; }
    void set_structure(std::string structure) {
        structure_ = structure;
        pairings.clear();
        pairs_struct.clear();
        unpaired_count_list.clear(); 
        is_pseudoknot_pair.clear();

        pairings.assign(structure_.size(), NULL_INDEX);
        is_pseudoknot_pair.assign(structure_.size(), false);
        pairs_struct.reserve(structure_.size());
        
        update_pairings();
        generate_unpaired_bases_count_list();
    }

    // [from, to)
    size_t get_unpaired_count(size_t from, size_t to) {
        if (from >= to) return 0;
        // todo: validate input
        return unpaired_count_list[to] - unpaired_count_list[from];
    }

   private:
    std::string name_;
    std::string sequence_;
    std::string structure_;
    std::vector<size_t>
        pairings;  // [4, -1, -1, 1] Number represents the index that base is paired to
    std::vector<Pair> pairs_struct;
    std::vector<size_t> unpaired_count_list;
    std::vector<bool> is_pseudoknot_pair;
    

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
        std::stack<size_t> brackets;
        std::stack<size_t> pseudoknots;
        size_t j;

        for (size_t i = 0; i < structure_.size(); i++) {
            switch (structure_[i]) {
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
                                                 sequence_ + "\nStructure: " + structure_);
                    }
                    j = brackets.top();
                    brackets.pop();
                    pairings[i] = j;
                    pairings[j] = i;
                    pairs_struct.emplace_back(i,j,false);
                    break;
                case ']':
                    if (pseudoknots.empty()) {
                        throw std::runtime_error("Structure in RNAEntry is invalid. \nSequence: " +
                                                 sequence_ + "\nStructure: " + structure_);
                    }
                    j = pseudoknots.top();
                    pseudoknots.pop();
                    pairings[i] = j;
                    pairings[j] = i;
                    is_pseudoknot_pair[i] = true;
                    is_pseudoknot_pair[j] = true;
                    pairs_struct.emplace_back(i,j,true);
                    break;
                default:
                    throw std::runtime_error(
                        "Character in RNAEntry's structure is invalid. \nInvalid Character: " +
                        std::string(1, structure_[i]) + "\nSequence: " + sequence_ +
                        "\nStructure: " + structure_);
            }
        }
        if (!brackets.empty() || !pseudoknots.empty()) {
            throw std::runtime_error("Unmatched opening brackets in RNA structure.\nSequence: " +
                                     sequence_ + "\nStructure: " + structure_);
        }
    }

    /**
     * @brief Creates a list indicating the number of unpaired bases up till that index
     */
    void generate_unpaired_bases_count_list() {
        size_t count = 0;
        size_t n = structure_.size();

        unpaired_count_list.assign(n + 1, 0);
        for (size_t i = 0; i < n; ++i) {
            count += (pairings[i] == NULL_INDEX);
            unpaired_count_list[i + 1] = count;
        }
    };
};
}  // namespace compute_energy