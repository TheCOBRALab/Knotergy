#pragma once

#include <iostream>
#include <stack>
#include <string>
#include <vector>

namespace ComputeEnergy {
struct RNAEntry {
    std::string name;
    std::string sequence;
    std::string structure;

    RNAEntry(std::string n, std::string s, std::string st)
        : name(std::move(n)), sequence(std::move(s)), structure(std::move(st)) {}

    // Default constructor (needed for vector resizing or default initialization)
    RNAEntry() = default;

    std::vector<int>& get_pairings() {
        // builtin_expect just tells the compiler that this case is unlikely to happen. It's used
        // for optimization
        if (__builtin_expect((structure != last_checked_structure), 0)) {
            update_pairings();
            last_checked_structure = structure;
        }
        return pairings;
    }

   private:
    std::vector<int>
        pairings;  // [4, -1, -1, 1] the number represents the index that base is paired to

    // Used to check if either structure or sequence has changed since last time the user got the
    // pairings
    std::string last_checked_structure = "";

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
        if (sequence.size() != structure.size()) {
            std::cerr << "Warning: Sequence and Structure are different sizes.\n"
                      << "Sequence: " << sequence << "\n"
                      << "Structure: " << structure << std::endl;
        }

        pairings.assign(sequence.size(), -1);  // pre-allocate pairings with -1
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
            if (!brackets.empty() || !pseudoknots.empty()) {
                throw std::runtime_error(
                    "Unmatched opening brackets in RNA structure.\nSequence: " + sequence +
                    "\nStructure: " + structure);
            }
        }
    }
};
}  // namespace ComputeEnergy