#pragma once

#include "preprocessing/ClosedRegion.hpp"
#include "preprocessing/RNAEntry.hpp"
#include "utils/common.hpp"

#include <uni_algo/ranges_grapheme.h>

#include <array>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>
namespace knotergy {

/**
 * @brief Immutable representation of an RNA entry with precomputed structural annotations.
 *
 * A ProcessedRNAEntry bundles:
 *  - The original RNA name, sequence, and structure string.
 *  - Base-pair index mappings (see RNAProcessor::compute_pair_table()).
 *  - Closed regions and their boundary pair_table (see RNAProcessor::compute_closed_regions()).
 *  - A prefix-sum array for fast unpaired-base counts (see RNAProcessor::compute_unpaired_counts())
 *  - Unpaired base counts within specified regions.
 *
 * This object is typically created by RNAProcessor::process_rna() and then queried by downstream
 * code.
 */
class ProcessedRNAEntry {
   private:
    std::string                   name_;                // RNA entry name.
    std::string                   raw_sequence_;        // Raw RNA input sequence
    std::string                   structure_;           // Dot-bracket RNA structure string.
    std::string                   sequence_;            // Unmodified RNA nucleotide sequence.
    std::vector<std::string_view> mod_sequence_views_;  // Modified RNA sequence.
    std::vector<size_t>           pair_table_;          // Base-pair indices for each position.
    std::vector<ClosedRegion>     closed_regions_;      // All closed regions in the structure.
    std::vector<size_t>           closed_regions_pair_table_;  // Closed regions boundary indicies.
    std::vector<int>              unpaired_prefix_sum_;  // Prefix-sum of unpaired-base counts.
    bool has_modified_bases_ = false;                    // Whether modified bases are present.

   public:
    /**
     * @brief Construct a ProcessedRNAEntry with precomputed annotations.
     *
     * @param name RNA entry name.
     * @param raw_sequence Raw RNA input sequence (may include modified bases).
     * @param structure Dot-bracket RNA structure string.
     * @param unmodified_sequence The unmodified RNA sequence (modified bases replaced).
     * @param mod_sequence_views The modified RNA sequence (raw sequence split into string_views per
     * base).
     * @param pair_table Base-pair indices for each position (NULL_INDEX for unpaired).
     * @param closed_regions List of closed regions detected in the structure.
     * @param closed_regions_pair_table Partner indices for closed-region boundaries.
     * @param unpaired_prefix_sum Prefix-sum array of unpaired-base counts (size = rna.size() + 1).
     * @param has_modified_bases Whether the RNA sequence contains modified bases.
     */
    ProcessedRNAEntry(std::string name, std::string raw_sequence, std::string structure,
                      std::string unmodified_sequence, std::vector<size_t> pair_table,
                      std::vector<ClosedRegion> closed_regions,
                      std::vector<size_t>       closed_regions_pair_table,
                      std::vector<int> unpaired_prefix_sum, bool has_modified_bases)
        : name_{std::move(name)},
          raw_sequence_{std::move(raw_sequence)},
          structure_{std::move(structure)},
          sequence_{std::move(unmodified_sequence)},
          pair_table_{std::move(pair_table)},
          closed_regions_{std::move(closed_regions)},
          closed_regions_pair_table_{std::move(closed_regions_pair_table)},
          unpaired_prefix_sum_{std::move(unpaired_prefix_sum)},
          has_modified_bases_{has_modified_bases} {
        if (has_modified_bases) {
            mod_sequence_views_ = compute_modified_sequence_views(raw_sequence_, structure_);
        }
    }

    /**
     * @brief Construct a ProcessedRNAEntry with precomputed annotations.
     *
     * This is a convenience constructor that extracts name, sequence, and structure
     * from an RNAEntry.
     *
     * @param rna Original RNAEntry (provides name, sequence, and structure).
     * @param unmodified_sequence The unmodified RNA sequence (modified bases replaced).
     * @param mod_sequence_views The modified RNA sequence (raw sequence split into string_views
     * @param pair_table Base-pair indices for each position (NULL_INDEX for unpaired).
     * @param closed_regions List of closed regions detected in the structure.
     * @param closed_regions_pair_table Partner indices for closed-region boundaries.
     * @param unpaired_prefix_sum Prefix-sum array of unpaired-base counts (size = rna.size() + 1).
     * @param has_modified_bases Whether the RNA sequence contains modified bases.
     */
    ProcessedRNAEntry(RNAEntry rna, std::string unmodified_sequence, std::vector<size_t> pair_table,
                      std::vector<ClosedRegion> closed_regions,
                      std::vector<size_t>       closed_regions_pair_table,
                      std::vector<int> unpaired_prefix_sum, bool has_modified_bases)
        : ProcessedRNAEntry(std::move(rna.name), std::move(rna.sequence), std::move(rna.structure),
                            std::move(unmodified_sequence), std::move(pair_table),
                            std::move(closed_regions), std::move(closed_regions_pair_table),
                            std::move(unpaired_prefix_sum), has_modified_bases) {}

    // Move constructor (ProcessedRNAEntry b = std::move(a);)
    ProcessedRNAEntry(ProcessedRNAEntry&& other) noexcept
        : name_(std::move(other.name_)),
          raw_sequence_(std::move(other.raw_sequence_)),
          structure_(std::move(other.structure_)),
          sequence_(std::move(other.sequence_)),
          mod_sequence_views_(),  // rebuild
          pair_table_(std::move(other.pair_table_)),
          closed_regions_(std::move(other.closed_regions_)),
          closed_regions_pair_table_(std::move(other.closed_regions_pair_table_)),
          unpaired_prefix_sum_(std::move(other.unpaired_prefix_sum_)),
          has_modified_bases_(other.has_modified_bases_) {
        if (has_modified_bases_) {
            mod_sequence_views_ = compute_modified_sequence_views(raw_sequence_, structure_);
        }
    }

    /**
     * @brief Move assignment operator (ProcessedRNAEntry a = std::move(b);)
     *
     * NOTE: We need to recompute mod_sequence_views_ if has_modified_bases_ is true, since
     * the moved reference string (raw_sequence_) may have been moved to a different memory
     * location.
     */
    ProcessedRNAEntry& operator=(ProcessedRNAEntry&& other) noexcept {
        if (this != &other) {
            name_                      = std::move(other.name_);
            raw_sequence_              = std::move(other.raw_sequence_);
            structure_                 = std::move(other.structure_);
            sequence_                  = std::move(other.sequence_);
            pair_table_                = std::move(other.pair_table_);
            closed_regions_            = std::move(other.closed_regions_);
            closed_regions_pair_table_ = std::move(other.closed_regions_pair_table_);
            unpaired_prefix_sum_       = std::move(other.unpaired_prefix_sum_);
            has_modified_bases_        = other.has_modified_bases_;
            if (has_modified_bases_) {
                mod_sequence_views_ = compute_modified_sequence_views(raw_sequence_, structure_);
            }
        }
        return *this;
    }

    /**
     * @brief Copy constructor (deep copy) (ProcessedRNAEntry a = b;)
     *
     * NOTE: We need to recompute mod_sequence_views_ if has_modified_bases_ is true, since the
     * modified sequence views are string_views into the raw_sequence string, and we need to ensure
     * they point to the correct memory location in the new copy.
     */
    ProcessedRNAEntry(const ProcessedRNAEntry& other)
        : name_(other.name_),
          raw_sequence_(other.raw_sequence_),
          structure_(other.structure_),
          sequence_(other.sequence_),
          pair_table_(other.pair_table_),
          closed_regions_(other.closed_regions_),
          closed_regions_pair_table_(other.closed_regions_pair_table_),
          unpaired_prefix_sum_(other.unpaired_prefix_sum_),
          has_modified_bases_(other.has_modified_bases_) {
        if (has_modified_bases_) {
            mod_sequence_views_ = compute_modified_sequence_views(raw_sequence_, structure_);
        }
    }

    /// @return The RNA entry's name.
    [[nodiscard]] const std::string& get_name() const { return name_; }

    /// @return The raw RNA input sequence. (includes modified bases)
    [[nodiscard]] const std::string& get_raw_sequence() const { return raw_sequence_; }

    /// @return The raw RNA structure (dot-bracket notation).
    [[nodiscard]] const std::string& get_structure() const { return structure_; }

    /// @return The unmodified RNA sequence (modified bases replaced with standard ones).
    [[nodiscard]] const std::string& get_sequence() const { return sequence_; }

    /// @return The modified RNA sequence (raw sequence split into string_views per base).
    [[nodiscard]] const std::vector<std::string_view>& get_modified_sequence() const {
        return mod_sequence_views_;
    }

    /// @return Base-pair indices for each position. See RNAProcessor::compute_pair_table for
    /// details
    [[nodiscard]] const std::vector<size_t>& get_pair_table() const { return pair_table_; }

    /// @return List of closed regions. See RNAProcessor::compute_closed_regions for details
    [[nodiscard]] const std::vector<ClosedRegion>& get_closed_regions() const {
        return closed_regions_;
    }

    /// Get closed-region indicies. See RNAProcessor::compute_cr_pair_table for details
    [[nodiscard]] const std::vector<size_t>& get_closed_regions_pair_table() const {
        return closed_regions_pair_table_;
    }

    /// @return Whether the RNA sequence contains modified bases or not
    [[nodiscard]] bool has_modified_bases() const { return has_modified_bases_; }

    /// @return The length of the RNA sequence/structure.
    [[nodiscard]] size_t size() const { return structure_.size(); }

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
    [[nodiscard]] int get_unpaired_count(size_t from, size_t to) const {
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
    [[nodiscard]] int get_unpaired_count(ClosedRegion closed_region) const {
        return get_unpaired_count(closed_region.begin, closed_region.end);
    }

    // grapheme is a user-perceived character, which may be multiple bytes. (e.g. Ψ is 2 bytes)
    // It's defined here since adding it to RNAProcessor would create a circular dependency
    [[nodiscard]] static std::vector<std::string_view> compute_modified_sequence_views(
        const std::string& sequence, const std::string& structure = "") {
        std::vector<std::string_view> out;
        out.reserve(sequence.size());

        // No structure provided = Use multi-byte parsing.
        // Sequence longer than structure = Use multi-byte parsing (Likely has multi-byte graphemes)
        // Multi-byte grapheme parsing is more expensive, so we want to avoid it if possible.
        if (structure.empty() || sequence.size() > structure.size()) {
            for (std::string_view g : una::views::grapheme::utf8(sequence)) {
                out.push_back(g);
            }
        } else {
            // Assumes each grapheme is a single character (Faster).
            for (size_t i = 0; i < sequence.size(); ++i) {
                unsigned char c = static_cast<unsigned char>(sequence[i]);

                // If the 0x80 bit is set, it's a multi-byte grapheme.
                // Fall back to multi-byte parsing for the entire sequence.
                if (c & 0x80) {
                    // No structure provided = Use multi-byte parsing.
                    return compute_modified_sequence_views(sequence);
                }
                out.emplace_back(sequence.data() + i, 1);
            }
        }

        return out;
    }
};

}  // namespace knotergy