#pragma once
#include <iostream>
#include <utility>

#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "LoopNode.hpp"
namespace knotergy {

/**
 * @brief Linked structure for tracking where the next paired base is.
 *
 * e.g.
 * ((..[[..))..]]
 * Linked list:
 * 0 -> 1 -> 4 -> 5 -> 8 -> 9 -> 12 -> 13 -> NULL_INDEX
 * 
 * Each node contains the position value, previous position, and next position.
 * This structure helps in navigating through base pairs when identifying bands.
 * 
 */
struct PairedBaseNode {
    size_t value = NULL_INDEX;  ///< Position value.
    size_t prev = NULL_INDEX;   ///< Previous position in the linked list.
    size_t next = NULL_INDEX;   ///< Next position in the linked list.
};

/**
 * @brief Identifies and constructs Band objects for pseudoknotted loops.
 *
 * This class analyzes RNA secondary structure pairing information to detect
 * pseudoknot bands. It identifies the four key positions (left_border, left_inner,
 * right_inner, right_border) that define each band.
 */
class BandFinder {
   public:
    BandFinder() = default;

    /**
     * @brief Find all bands within a specified region of an RNA structure.
     *
     * @param left_bound Left boundary of the region to search.
     * @param right_bound Right boundary of the region to search.
     * @param loop_type The type of loop being analyzed.
     * @param pairings Base-pair index mapping for the structure.
     * @param cr_pairings Closed region pairing indices.
     * @return Vector of Band objects found in the region.
     */
    static std::vector<Band> find_bands(const size_t& left_bound, const size_t& right_bound,
                                        const LoopType& loop_type,
                                        const std::vector<size_t>& pairings,
                                        const std::vector<size_t>& cr_pairings, bool simple_bands = false);

    /**
     * @brief Find all bands for a given loop node.
     *
     * Convenience method that extracts bounds from the loop node.
     *
     * @param node The loop node to analyze.
     * @param processed_rna The processed RNA entry with pairing information.
     * @return Vector of Band objects found for this loop node.
     */
    static std::vector<Band> find_bands(const LoopNode& node,
                                        const ProcessedRNAEntry& processed_rna);

   private:
    /**
     * @brief Extend a stem to find the inner band positions.
     *
     * Follows consecutive stacking pairs from the initial positions to find
     * the innermost positions of a band's left and right arms.
     *
     * @param i_prime Left inner position (modified in place).
     * @param j_prime Right inner position (modified in place).
     * @param aux_bands Map of band links for navigation.
     * @param pairings Base-pair index mapping.
     * @return True if the stem was successfully extended.
     */
    static bool extend_stem(size_t& i_prime, size_t& j_prime,
                            const std::unordered_map<size_t, PairedBaseNode>& aux_bands,
                            const std::vector<size_t>& pairings);

    /**
     * @brief Generate a linked structure of potential band boundaries.
     *
     * Creates a map linking positions that could form band boundaries.
     *
     * @param left_bound Left boundary of the region.
     * @param right_bound Right boundary of the region.
     * @param pairings Base-pair index mapping.
     * @param cr_pairings Closed region pairing indices.
     * @return Map of band links indexed by position.
     */
    static std::unordered_map<size_t, PairedBaseNode> const generate_paired_base_links(
        const size_t& left_bound, const size_t& right_bound, const std::vector<size_t>& pairings,
        const std::vector<size_t>& cr_pairings);

    /**
     * @brief Generate band links for a given loop node.
     *
     * Convenience method that extracts bounds from the loop node.
     *
     * @param node The loop node to analyze.
     * @param processed_entry The processed RNA entry with pairing information.
     * @return Map of band links indexed by position.
     */
    static std::unordered_map<size_t, PairedBaseNode> const generate_paired_base_links(
        const LoopNode& node, const ProcessedRNAEntry& processed_entry);
};
}  // namespace knotergy