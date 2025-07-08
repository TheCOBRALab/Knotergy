#pragma once

#include <array>
#include <iostream>
#include <stack>
#include <string>
#include <vector>

#include "../helpers/common.hpp"
#include "ClosedRegion.hpp"
#include "RNAEntry.hpp"
namespace knotergy {

class RNAProcessedEntry {
   public:
    // Constructor for RNAEntry with a name
    explicit RNAProcessedEntry(const RNAEntry& entry);
    explicit RNAProcessedEntry(std::string name, std::string sequence, std::string structure);
    explicit RNAProcessedEntry(std::string sequence, std::string structure);

    RNAProcessedEntry() = default;

    [[nodiscard]] const std::string& get_name() const;
    [[nodiscard]] const std::string& get_sequence() const;
    [[nodiscard]] const std::string& get_structure() const;
    [[nodiscard]] const std::vector<size_t>& get_pairings() const;
    [[nodiscard]] const std::vector<ClosedRegion>& get_closed_regions() const;
    [[nodiscard]] size_t size() const;

    // [from, to)
    [[nodiscard]] int get_unpaired_count(size_t from, size_t to) const;
    [[nodiscard]] int get_unpaired_count(ClosedRegion cr) const;

   private:
    // pairings represents the indicies where base is paired to. e.g. (..) = [3, -1, -1, 0]
    std::vector<size_t> pairings_;

    // All closed and weakly closed regions
    // Closed region = a region where every base that is opened is also closed
    // e.g. ((..))..([..)]. -> (0, 5), (1, 4), (8, 13)
    std::vector<ClosedRegion> closed_regions_;

    // Computes all the unpaired count as a prefix sum with an added 0 at the start
    // e.g. .(.). -> 0, 1, 1, 2, 2, 3
    // unpaired_prefix_sum_[j] - unpaired_prefix_sum_[i] = unpaired count in range [i,j)
    std::vector<int> unpaired_prefix_sum_;

    // RNA name, sequence and structure
    RNAEntry rna_;

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
    std::vector<size_t> compute_pairings();

    std::vector<ClosedRegion> compute_closed_regions();

    /**
     * @brief Creates a list indicating the number of unpaired bases up till that index
     */
    std::vector<int> compute_unpaired_counts();
};

}  // namespace knotergy