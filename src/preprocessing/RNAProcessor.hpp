#pragma once

#include "io/output/colors.hpp"
#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "preprocessing/ClosedRegion.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
#include "preprocessing/RNAEntry.hpp"

#include <array>
#include <vector>
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
 * - Unpaired positions are marked with NULL_INDEX (defined as static_cast<std::size_t>(-1)).
 *   This evaluates to the maximum std::size_t value; do not assume it is negative.
 * - Supported bracket types: (), [], {}, <>.
 */
class RNAProcessor {
   public:
    RNAProcessor();

    /// High-level pipeline that produces a processed entry from raw RNA input.
    static ProcessedRNAEntry process_rna(RNAEntry rna, const all_mod_params& modified_params = {});

    /**
     * @brief Compute base-pair indices from a structure string.
     *
     * Scans an RNA structure (dot-bracket notation) and returns, for each position i,
     * the index j of its paired base, or NULL_INDEX if unpaired.
     *
     * Examples:
     * - "(..)" → pair_table = [3, NULL_INDEX, NULL_INDEX, 0]
     * - "([..)]" → pair_table = [4, 5, NULL_INDEX, NULL_INDEX, 0, 1]
     *
     * Notes:
     * - Indices are 0-based.
     * - Supported bracket pairs: (), [], {}, <>. Different bracket types may be nested or
     * intermixed.
     * - NULL_INDEX (defined as static_cast<std::size_t>(-1)) marks unpaired bases.
     *   This evaluates to the maximum std::size_t value because std::size_t is unsigned.
     *
     * @param structure Dot-bracket RNA structure string.
     * @param unmodified_sequence (Optional) The unmodified RNA sequence corresponding to the
     * structure.
     * @param mod_sequence (Optional) Modified RNA sequence (raw sequence split into string_views
     * per base).
     * @return std::vector<std::size_t> of length rna.size(), where pair_table[i] is the index of
     * i's partner, or NULL_INDEX if i is unpaired.
     *
     * @throws std::runtime_error If the structure is malformed (e.g., unbalanced/mismatched
     * brackets).
     * @warning Invalid base *types* (e.g., A–A) are reported via warnings but do not throw.
     */
    [[nodiscard]] static std::vector<std::size_t> compute_pair_table(
        const std::string& structure, const std::string& unmodified_sequence = "",
        const std::vector<std::string_view>& mod_sequence = {});

    /**
     * @brief Compute base-pair indices from a structure string.
     *
     * Convenience overload: extracts the structure from an RNAEntry and passes it to the main
     * compute_pair_table() method.
     *
     * @param rna RNAEntry containing at least the structure
     * @param unmodified_sequence (Optional) The unmodified RNA sequence corresponding to the
     *structure.
     * @param mod_sequence (Optional) Modified RNA sequence (raw sequence split into string_views
     *per base).
     **/
    [[nodiscard]] static std::vector<std::size_t> compute_pair_table(
        const RNAEntry& rna, const std::string& unmodified_sequence = "",
        const std::vector<std::string_view>& mod_sequence = {});

    /**
     * @brief Identify all closed regions in the structure.
     *
     * Parses the base-pair vector returned by compute_pair_table() and returns every closed region.
     * See ClosedRegion.hpp for the definition and semantics of a closed region.
     *
     * Output is sorted by the start index of each closed region
     *
     * Example:
     *  pair_table = [3, NULL_INDEX, NULL_INDEX, 0] → closed_regions = [ClosedRegion(0, 3)]
     *
     * @param pair_table Base-pair indices as returned by compute_pair_table().
     * @param number_of_pairs (Optional) If known, the number of base pairs can be passed to
     * preallocate the result vector.
     * @return std::vector<ClosedRegion> containing all detected closed regions.
     */
    [[nodiscard]] static std::vector<ClosedRegion> compute_closed_regions(
        const std::vector<std::size_t>& pair_table, std::size_t number_of_pairs = 32);

    /**
     * @brief Prefix-sum of unpaired-base counts.
     *
     * Computes a vector U of length (rna_size + 1) where:
     *   - U[0] = 0
     *   - For k in [1, rna_size], U[k] is the number of unpaired bases in pair_table[0..k-1].
     *
     * Example:
     *   pair_table = [3, NULL_INDEX, NULL_INDEX, 0]
     *   → U = [0, 0, 1, 2, 2]
     *
     * To count unpaired bases in a half-open interval [start, end):
     *   unpaired = U[end] - U[start]
     *
     * @param pair_table Base-pair indices from RNAProcessor::compute_pair_table().
     * @return std::vector<int> of size rna_size + 1 with cumulative unpaired counts.
     */
    [[nodiscard]] static std::vector<int> compute_unpaired_counts(
        const std::vector<std::size_t>& pair_table);

    /**
     * @brief Compute the unmodified RNA sequence from the modified sequence and parameters.
     *
     * This is necessary as the energy is first computed on the unmodified sequence,
     * and then adjusted based on the modified bases.
     *
     * Example:
     * modified_sequence = ["6", "U", "G", "C", "P"] -> unmodified_sequence = "AUGCU"
     */
    [[nodiscard]] static std::string compute_unmodified_sequence(
        const std::vector<std::string_view>& modified_sequence_views, const all_mod_params& params,
        const std::size_t rna_length, bool& has_modified_bases);

    /**
     * @brief Check if a base is an unmodified base (A, U, G, C, T, N).
     *
     * This is used to validate the input and ensure that modified bases are properly handled.
     *
     * Example:
     * is_unmodified_base('A') -> true
     * is_unmodified_base('6') -> false
     *
     * @param base The base to check, as a string_view or char.
     * @return true if the base is unmodified, false otherwise.
     */
    [[nodiscard]] static bool is_unmodified_base(const std::string_view& base);
    [[nodiscard]] static bool is_unmodified_base(char base);

   private:
    // Lookup table for unmodified bases
    [[nodiscard]] static constexpr bool is_unmodified_base(unsigned char c) noexcept {
        switch (c) {
            case 'A':
            case 'C':
            case 'G':
            case 'N':
            case 'T':
            case 'U': return true;
            default:  return false;
        }
    }

    // Function to check if two bases can pair according to RNA base-pairing rules.
    [[nodiscard]] static constexpr bool can_pair(char left, char right) {
        switch (left) {
            case 'A': return right == 'U' || right == 'T';
            case 'U': return right == 'A' || right == 'G';
            case 'G': return right == 'C' || right == 'U' || right == 'T';
            case 'C': return right == 'G';
            case 'T': return right == 'A' || right == 'G';
            case 'N': return false;  // Should not pair with anything
            default:  return false;
        }
    };
};

}  // namespace knotergy