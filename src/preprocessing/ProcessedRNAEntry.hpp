#pragma once

#include <array>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "../io/common.hpp"
#include "ClosedRegion.hpp"
#include "RNAEntry.hpp"
namespace knotergy {


/**
 * @brief Immutable representation of an RNA entry with precomputed structural annotations.
 *
 * A ProcessedRNAEntry bundles:
 *  - The original RNA name, sequence, and structure string.
 *  - Base-pair index mappings (see RNAProcessor::compute_pairings()).
 *  - Closed regions and their boundary pairings (see RNAProcessor::compute_closed_regions()).
 *  - A prefix-sum array for fast unpaired-base counts (see RNAProcessor::compute_unpaired_counts()).
 *  - Unpaired base counts within specified regions.
 *
 * This object is typically created by RNAProcessor::process_rna() and then queried by downstream code.
 */
class ProcessedRNAEntry {
   public:
   /**
     * @brief Construct a ProcessedRNAEntry with precomputed annotations.
     *
     * @param rna Original RNA entry (provides name, sequence, and structure).
     * @param pairings Base-pair indices for each position (NULL_INDEX for unpaired).
     * @param closed_regions List of closed regions detected in the structure.
     * @param closed_regions_pairings Partner indices for closed-region boundaries.
     * @param unpaired_prefix_sum Prefix-sum array of unpaired-base counts (size = rna.size() + 1).
     */
        explicit ProcessedRNAEntry(
            RNAEntry rna,
            std::string unmodified_sequence,
            std::vector<std::string_view> mod_sequence_views,
            std::vector<size_t> pairings,
            std::vector<ClosedRegion> closed_regions,
            std::vector<size_t> closed_regions_pairings,
            std::vector<int> unpaired_prefix_sum,
            bool has_modified_bases)
            : name_{rna.name},
            sequence_{unmodified_sequence},
            mod_sequence_views_{mod_sequence_views},
            raw_sequence_{std::move(rna.sequence)},
            structure_{rna.structure},
            pairings_{pairings},
            closed_regions_{closed_regions},
            closed_regions_pairings_{closed_regions_pairings},
            unpaired_prefix_sum_{unpaired_prefix_sum},
            has_modified_bases_{has_modified_bases} {}


    /// @return The RNA entry's name.
    const std::string& get_name() const { return name_; }

    /// @return The raw RNA input sequence. (includes modified bases)
    const std::string& get_raw_sequence() const { return raw_sequence_; }
    
    /// @return The unmodified RNA sequence (modified bases replaced with standard ones).
    const std::string& get_sequence() const { return sequence_; }

    /// @return The modified RNA sequence (raw sequence split into string_views per base).
    const std::vector<std::string_view>& get_modified_sequence() const { return mod_sequence_views_; }

    /// @return The raw RNA structure (dot-bracket notation).
    const std::string& get_structure() const { return structure_; }

    /// Get the base-pair index vector. See RNAProcessor::compute_pairings() for details.
    const std::vector<size_t>& get_pairings() const { return pairings_; }

    ///Get list of closed regions. See RNAProcessor::compute_closed_regions() for details.
    const std::vector<ClosedRegion>& get_closed_regions() const { return closed_regions_; }

    /// Get closed-region indicies. See RNAProcessor::compute_cr_pairings() for details.
    const std::vector<size_t>& get_closed_regions_pairings() const {
        return closed_regions_pairings_;
    }

    /// @return Whether the RNA sequence contains modified bases.
    bool has_modified_bases() const { return has_modified_bases_; }

    /// @return The length of the RNA sequence/structure.
    size_t size() const { return structure_.size(); }

    /**
     * @brief Compute the number of unpaired bases in a half-open interval [from, to).
     *
     * Uses the prefix-sum vector computed at construction.
     * 
     * @see RNAProcessor::compute_unpaired_counts() for more details.
     *
     * @param from Start index (inclusive).
     * @param to End index (exclusive).
     * @return Number of unpaired bases between from and to.
     * @throws std::out_of_range if indices are out of bounds.
     */
    int get_unpaired_count(size_t from, size_t to) const {
        if (from >= unpaired_prefix_sum_.size() || to > unpaired_prefix_sum_.size()) {
            throw std::out_of_range("Index out of range in get_unpaired_count");
        }
        if (from >= to) return 0;

        return unpaired_prefix_sum_[to] - unpaired_prefix_sum_[from];
    }

    /**
     * @brief Compute the number of unpaired bases inside a closed region.
     *
     * Equivalent to get_unpaired_count(region.begin, region.end).
     */
    int get_unpaired_count(ClosedRegion closed_region) const {
        return get_unpaired_count(closed_region.begin, closed_region.end);
    }

   private:
    const std::string name_;                                  // RNA entry name.
    const std::string sequence_;                              // Unmodified RNA nucleotide sequence.
    const std::vector<std::string_view> mod_sequence_views_;  // Modified RNA sequence.
    std::string raw_sequence_;                                // Raw RNA input sequence. (not const to allow move) (move required for string_view vector)
    const std::string structure_;                             // Dot-bracket RNA structure string.
    const std::vector<size_t> pairings_;                      // Base-pair indices for each position.
    const std::vector<ClosedRegion> closed_regions_;          // All closed regions in the structure.
    const std::vector<size_t> closed_regions_pairings_;       // Closed regions boundary indicies.
    const std::vector<int> unpaired_prefix_sum_;              // Prefix-sum of unpaired-base counts.
    const bool has_modified_bases_ = false;                   // Whether modified bases are present.
};

}  // namespace knotergy