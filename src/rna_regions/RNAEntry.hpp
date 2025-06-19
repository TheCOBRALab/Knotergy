#pragma once

#include <iostream>
#include <stack>
#include <string>
#include <vector>

#include "../helpers/common.hpp"

namespace compute_energy {

struct Region {
    size_t begin{};
    size_t end{};
    bool pseudoknotted = false;

    Region() = default;
    Region(size_t b, size_t e) : begin(b), end(e) {}
    Region(size_t b, size_t e, bool p) : begin(b), end(e), pseudoknotted(p) {}

    bool operator==(const Region& rhs) const { return begin == rhs.begin && end == rhs.end; }
};

// Lets you print out the Region (overloading the << operator )
inline std::ostream& operator<<(std::ostream& os, const Region& region) {
    os << "Region(" << region.begin << ", " << region.end << ") pk: " << region.pseudoknotted;
    return os;
}

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
    [[nodiscard]] const std::vector<Region>& get_closed_regions() const { return closed_regions; };
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
        unpaired_count_list.clear();
        pairings.assign(structure_.size(), NULL_INDEX);

        update_pairings();
        generate_unpaired_bases_count_list();
        compute_closed_regions();
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
    // pairings represents the indicies where base is paired to. e.g. (..) = [3, -1, -1, 0]
    std::vector<size_t> pairings;
    std::vector<Region> closed_regions;
    std::vector<size_t> unpaired_count_list;

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
    }

    void compute_closed_regions() {
        closed_regions.clear();
        std::stack<Region> stack;
        const size_t n = pairings.size();

        for (size_t i = 0; i < n; ++i) {
            size_t bp = pairings[i];
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

// Operator overloading to output all closed regions (cout << Region)
inline std::ostream& operator<<(std::ostream& os, const RNAEntry& entry) {
    os << "Closed Regions:\n";
    for (const Region& r : entry.get_closed_regions()) {
        os << "Region from " << r.begin << " to " << r.end << ": ";
        os << entry.get_sequence().substr(r.begin, r.end - r.begin + 1) << std::endl;
    }
    return os;
}
}  // namespace compute_energy