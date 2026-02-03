#pragma once

#include <array>
#include <vector>

#include "../io/ViennaParams.hpp"
#include "ClosedRegion.hpp"
#include "ProcessedRNAEntry.hpp"
#include "RNAEntry.hpp"
namespace knotergy {

/**
 * @brief Utilities for parsing RNA secondary-structure strings.
 *
 * This class provides helpers to:
 * - compute base-pair indices from a dot-bracket–style structure,
 * - extract closed regions,
 * - and precompute convenience arrays (e.g., unpaired-base counts).
 *
 * Conventions:
 * - Indices are 0-based.
 * - Unpaired positions are marked with NULL_INDEX (defined as static_cast<size_t>(-1)).
 *   This evaluates to the maximum size_t value; do not assume it is negative.
 * - Supported bracket types: (), [], {}, <>.
 */
class RNAProcessor {
   public:
    RNAProcessor();

    /// High-level pipeline that produces a processed entry from raw RNA input.
    static ProcessedRNAEntry process_rna(
        RNAEntry rna, const std::vector<modified_base_params>& modified_params = {});

    /**
     * @brief Compute base-pair indices from a structure string.
     *
     * Scans an RNA structure (dot-bracket notation) and returns, for each position i,
     * the index j of its paired base, or NULL_INDEX if unpaired.
     *
     * Examples:
     * - "(..)" → pairings = [3, NULL_INDEX, NULL_INDEX, 0]
     * - "([..)]" → pairings = [4, 5, NULL_INDEX, NULL_INDEX, 0, 1]
     *
     * Notes:
     * - Indices are 0-based.
     * - Supported bracket pairs: (), [], {}, <>. Different bracket types may be nested or
     * intermixed.
     * - NULL_INDEX (defined as static_cast<size_t>(-1)) marks unpaired bases.
     *   This evaluates to the maximum size_t value because size_t is unsigned.
     *
     * @param rna RNAEntry containing the structure to parse.
     * @param unmodified_sequence The unmodified RNA sequence corresponding to the structure.
     * @return std::vector<size_t> of length rna.size(), where pairings[i] is the index of i's
     * partner, or NULL_INDEX if i is unpaired.
     *
     * @throws std::runtime_error If the structure is malformed (e.g., unbalanced/mismatched
     * brackets).
     * @warning Invalid base *types* (e.g., A–A) are reported via warnings but do not throw.
     */
    [[nodiscard]] static std::vector<size_t> compute_pairings(
        const RNAEntry& rna, const std::string& unmodified_sequence,
        const std::vector<std::string_view>& mod_sequence = {});

    /**
     * @brief Identify all closed regions in the structure.
     *
     * Parses the base-pair vector returned by compute_pairings() and returns every closed region.
     * See ClosedRegion.hpp for the definition and semantics of a closed region.
     *
     * @param pairings Base-pair indices as returned by compute_pairings().
     * @return std::vector<ClosedRegion> containing all detected closed regions.
     */
    [[nodiscard]] static std::vector<ClosedRegion> compute_closed_regions(
        const std::vector<size_t>& pairings);

    /**
     * @brief Build a partner-index vector for closed-region boundaries.
     *
     * Similar to RNAProcessor::compute_pairings(), but stores only the closed regions,
     * not every individual base pair.
     *
     * Example:
     *   [ClosedRegion(0, 5), ClosedRegion(2, 4)] → [5, NULL_INDEX, 4, NULL_INDEX, 2, 0]
     *
     * This is useful for quickly skipping over already-processed closed regions.
     *
     * @param closed_regions All closed regions in the structure.
     * @param rna_size The total length of the structure (used to preallocate the result).
     * @return std::vector<size_t> of length rna_size where entries are boundary partners or
     * NULL_INDEX.
     *
     * @throws std::runtime_error If any closed-region index exceeds rna_size.
     */
    [[nodiscard]] static std::vector<size_t> compute_cr_pairings(
        const std::vector<ClosedRegion>& closed_regions, const size_t& rna_size);

    /**
     * @brief Prefix-sum of unpaired-base counts.
     *
     * Computes a vector U of length (rna_size + 1) where:
     *   - U[0] = 0
     *   - For k in [1, rna_size], U[k] is the number of unpaired bases in pairings[0..k-1].
     *
     * Example:
     *   pairings = [3, NULL_INDEX, NULL_INDEX, 0]
     *   → U = [0, 0, 1, 2, 2]
     *
     * To count unpaired bases in a half-open interval [start, end):
     *   unpaired = U[end] - U[start]
     *
     * @param pairings Base-pair indices from RNAProcessor::compute_pairings().
     * @return std::vector<int> of size rna_size + 1 with cumulative unpaired counts.
     */
    [[nodiscard]] static std::vector<int> compute_unpaired_counts(
        const std::vector<size_t>& pairings);

    [[nodiscard]] static std::string compute_unmodified_sequence(
        const std::vector<std::string_view>& modified_sequence_views,
        const std::vector<modified_base_params>& params, const size_t rna_length,
        bool& has_modified_bases);

    [[nodiscard]] static bool is_unmod_base(const std::string_view& base);

   private:
    // Lookup table for unmodified bases
    static constexpr std::array<uint8_t, 256> unmod_lookup = [] {
        std::array<uint8_t, 256> t{};
        t[static_cast<unsigned char>('A')] = 1;
        t[static_cast<unsigned char>('U')] = 1;
        t[static_cast<unsigned char>('G')] = 1;
        t[static_cast<unsigned char>('C')] = 1;
        t[static_cast<unsigned char>('T')] = 1;
        t[static_cast<unsigned char>('N')] = 1;
        return t;
    }();
};

}  // namespace knotergy