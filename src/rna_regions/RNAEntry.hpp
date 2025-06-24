#pragma once

#include <array>
#include <iostream>
#include <stack>
#include <string>
#include <vector>

#include "../helpers/common.hpp"
#include "ClosedRegion.hpp"

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

    // ----------------------------------------- Getters -----------------------------------------
    [[nodiscard]] const std::string& get_name() const { return name_; }
    [[nodiscard]] const std::string& get_sequence() const { return sequence_; }
    [[nodiscard]] const std::string& get_structure() const { return structure_; }
    [[nodiscard]] const std::vector<ClosedRegion>& get_closed_regions() const {
        return closed_regions_;
    };
    [[nodiscard]] const std::vector<size_t>& get_pairings() const {
        if (sequence_.size() != structure_.size()) {
            std::cerr << "Warning: Sequence and Structure are different sizes.\n"
                      << "Sequence: " << sequence_ << "\n"
                      << "Structure: " << structure_ << std::endl;
        }
        return pairings_;
    }
    [[nodiscard]] const std::vector<std::array<size_t, 4>>& get_bands() const { return bands_; }

    // ---------------------------------------- Setters -----------------------------------------
    void set_name(std::string name) { name_ = name; }
    void set_sequence(std::string sequence) { sequence_ = sequence; }
    void set_structure(std::string structure) {
        structure_ = structure;

        pairings_ = update_pairings();
        unpaired_count_list_ = generate_unpaired_bases_count_list();
        closed_regions_ = compute_closed_regions();
        bands_ = find_bands();
    }

    // [from, to)
    int get_unpaired_count(size_t from, size_t to) const {

        if (from >= unpaired_count_list_.size() || to > unpaired_count_list_.size()) {
            throw std::out_of_range("Index out of range in get_unpaired_count");
        }

        if (from >= to)
            return 0;
        
        // Return the difference in unpaired counts between the two indices
        return unpaired_count_list_[to] - unpaired_count_list_[from];
    }

    int get_unpaired_count(ClosedRegion cr) const {
        return get_unpaired_count(cr.begin, cr.end);
    }

    // --------------------------------- Proccess Structure ---------------------------------
   private:
    std::string name_;
    std::string sequence_;
    std::string structure_;

    // pairings represents the indicies where base is paired to. e.g. (..) = [3, -1, -1, 0]
    std::vector<size_t> pairings_;

    std::vector<ClosedRegion> closed_regions_;
    std::vector<std::array<size_t, 4>> bands_;
    std::vector<int> unpaired_count_list_;

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
    std::vector<size_t> update_pairings() {
        std::stack<size_t> brackets;
        std::stack<size_t> pseudoknots;
        size_t j;
        std::vector<size_t> pairings(structure_.size(), NULL_INDEX);

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
        return pairings;
    }

    std::vector<ClosedRegion> compute_closed_regions() {
        std::vector<ClosedRegion> closed_regions;
        std::stack<ClosedRegion> stack;
        const size_t n = pairings_.size();

        for (size_t i = 0; i < n; ++i) {
            size_t bp = pairings_[i];
            if (bp == NULL_INDEX) continue;  // unpaired

            // ───── OPENING BASE: i < bp ────────────────────────────
            if (i < bp) {
                stack.push({i, bp});
                continue;
            }

            // ───── CLOSING BASE: bp < i ────────────────────────────
            size_t largest_right = i;  // rightmost boundary seen

            // if crossing (pseudoknotted), find right end of closed region
            while (!stack.empty() && stack.top().begin > bp) {
                stack.top().pseudoknotted = true;
                largest_right = std::max(largest_right, stack.top().end);
                stack.pop();
            }

            if (stack.empty()) continue;  // if unbalanced (should never happen)

            // extend region if needed
            stack.top().end = std::max(largest_right, stack.top().end);

            // region finished?
            if (i == stack.top().end) {
                closed_regions.push_back(stack.top());
                stack.pop();
            }
        }
        return closed_regions;
    }

    /**
     * @brief Creates a list indicating the number of unpaired bases up till that index
     */
    std::vector<int> generate_unpaired_bases_count_list() {
        int count = 0;
        size_t n = structure_.size();
        std::vector<int> unpaired_count_list;

        unpaired_count_list.assign(n + 1, 0);
        for (size_t i = 0; i < n; ++i) {
            count += (pairings_[i] == NULL_INDEX);
            unpaired_count_list[i + 1] = count;
        }
        return unpaired_count_list;
    };

    // ─────────────────────────────────────────────────────────────
    //  Helper:  extend one step along a perfectly stacked stem.
    //  Return true while i+1 pairs j-1.
    bool extend_stem(size_t& il, size_t& jr) const {
        const size_t n = pairings_.size();

        size_t ip = il + 1;
        size_t jp = jr - 1;

        // find next valid base pair
        while (ip < n && (pairings_[ip] == NULL_INDEX)) {
            ++ip;
        }

        // find next valid base pair
        while (jp > 0 && (pairings_[jp] == NULL_INDEX)) {
            --jp;
        }

        if (ip < jp && pairings_[ip] == jp) {
            il = ip;
            jr = jp;
            return true;
        }
        return false;
    }

    // ─────────────────────────────────────────────────────────────
    //  Main: find all bands inside the whole molecule
    //  Result: vector<array<size_t,4>>  storing {i, i′, j′, j}
    std::vector<std::array<size_t, 4>> find_bands() const {
        const size_t n = pairings_.size();

        // pair is already used in a previous band?
        std::vector<bool> done(n, false);
        std::vector<std::array<size_t, 4>> bands;

        for (size_t i = 0; i < n; ++i) {
            // left border of an unexplored base pair?
            if (done[i] || pairings_[i] == NULL_INDEX || pairings_[i] < i) continue;

            size_t j = pairings_[i];
            size_t il = i;  // will become i′
            size_t jr = j;  // will become j′

            // walks the stem until last base pair in the band is found
            while (extend_stem(il, jr)) {
            }

            // record the band
            bands.push_back(std::array<size_t, 4>{i, il, jr, j});

            // mark every base that is now part of a finished band as done
            for (size_t k = i; k <= il; ++k) {
                done[k] = true;
                if (pairings_[k] != NULL_INDEX) done[pairings_[k]] = true;
            }

            for (size_t k = jr; k <= j; ++k) {
                done[k] = true;
                if (pairings_[k] != NULL_INDEX) done[pairings_[k]] = true;
            }

            // fast-forward outer loop (go to next base pair)
            i = il;
        }
        return bands;
    }
};

// Operator overloading to output all closed regions (cout << Region)
inline std::ostream& operator<<(std::ostream& os, const RNAEntry& entry) {
    os << "Closed Regions:\n";
    for (const ClosedRegion& r : entry.get_closed_regions()) {
        os << "Region from " << r.begin << " to " << r.end << ": ";
        os << entry.get_sequence().substr(r.begin, r.end - r.begin + 1) << std::endl;
    }
    return os;
}
}  // namespace compute_energy