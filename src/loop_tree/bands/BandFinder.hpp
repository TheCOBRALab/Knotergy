#pragma once
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

#include <iostream>
#include <utility>
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
    size_t prev = NULL_INDEX;  ///< Previous position in the linked list.
    size_t next = NULL_INDEX;  ///< Next position in the linked list.
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
     * @param pair_table Base-pair index mapping for the structure.
     * @return Vector of Band objects found in the region.
     */
    [[nodiscard]] static std::vector<Band> find_bands(size_t cr_start, size_t cr_end,
                                                      LoopType loop_type,
                                                      std::vector<PairedBaseNode>& aux_bands,
                                                      const std::vector<size_t>& pair_table,
                                                      const std::vector<LoopNode*>& node_table);

    /**
     * @brief Find all bands for a given loop node.
     *
     * Convenience method that extracts bounds from the loop node.
     *
     * @param node The loop node to analyze.
     * @param processed_rna The processed RNA entry with pairing information.
     * @return Vector of Band objects found for this loop node.
     */
    [[nodiscard]] static std::vector<Band> find_bands(const LoopNode& node,
                                                      std::vector<PairedBaseNode>& aux_bands,
                                                      const ProcessedRNAEntry& processed_rna,
                                                      const std::vector<LoopNode*>& node_table);

   private:
    /**
     * @brief Extend a stem to find the inner band positions.
     *
     * Follows consecutive stacking pairs from the initial positions to find
     * the innermost positions of a band's left and right arms.
     *
     * @param i Left border of band.
     * @param j Right border of band.
     * @param aux_bands Map of band links for navigation.
     * @param pair_table Base-pair index mapping.
     * @return A pair containing the extended inner positions.
     */
    [[nodiscard]] static std::pair<size_t, size_t> find_stem_inner_indices(
        size_t band_start, size_t band_end, const std::vector<PairedBaseNode>& aux_bands,
        const std::vector<size_t>& pair_table);

    /**
     * @brief Generate a linked structure of potential band boundaries.
     *
     * Links consecutive paired bases within a closed region.
     * Nested children are skipped and are processed separately.
     * Each closed region is treated as an isolated unit,
     * so links are only generated within the bounds of the closed region.
     *
     * However this function only proesses one closed region at a time.
     * The example below is how the aux_bands would look after processing the entire structure.
     *
     * e.g.
     * Full Structure: ([.(.).([)].)]
     *
     * idx
     * 0 : prev = NULL_INDEX, next = 1
     * 1 : prev = 0,          next = 12
     * 2 : prev = NULL_INDEX, next = NULL_INDEX (unpaired)
     * 3 : prev = NULL_INDEX, next = NULL_INDEX (pseudoknot-free)
     * 4 : prev = NULL_INDEX, next = NULL_INDEX (unpaired)
     * 5 : prev = NULL_INDEX, next = NULL_INDEX (pseudoknot-free)
     * 6 : prev = NULL_INDEX, next = NULL_INDEX (unpaired)
     * 7 : prev = NULL_INDEX, next = 8
     * 8 : prev = 7,          next = 9
     * 9 : prev = 8,          next = 10
     * 10: prev = 9,          next = NULL_INDEX
     * 11: prev = NULL_INDEX, next = NULL_INDEX (unpaired)
     * 12: prev = 1,          next = 13
     * 13: prev = 12,         next = NULL_INDEX
     *
     *
     * @param left_bound Left boundary of the region.
     * @param right_bound Right boundary of the region.
     * @param pair_table Base-pair index mapping.
     * @param node_table Vector of all loop nodes.
     * @return Vector of PairedBaseNode objects.
     */
    static void generate_paired_base_links(size_t cr_start, size_t cr_end,
                                           std::vector<PairedBaseNode>& aux_bands,
                                           const std::vector<size_t>& pair_table,
                                           const std::vector<LoopNode*>& node_table);

    /**
     * @brief Generate band links for a given loop node.
     *
     * Convenience method that extracts bounds from the loop node.
     *
     * @param node The loop node to analyze.
     * @param processed_entry The processed RNA entry with pairing information.
     * @return Vector of PairedBaseNode objects.
     */
    static void generate_paired_base_links(const LoopNode& node,
                                           std::vector<PairedBaseNode>& aux_bands,
                                           const ProcessedRNAEntry& processed_entry,
                                           const std::vector<LoopNode*>& node_table);
};
}  // namespace knotergy