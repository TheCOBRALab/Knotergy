#pragma once
#include "loop_tree/LoopNode.hpp"
#include "loop_tree/bands/Band.hpp"
#include "preprocessing/ClosedRegion.hpp"

#include <vector>

namespace knotergy {

class BandBuilder {
   public:
    /**
     * @brief Construct a Band from boundary positions and pairing information.
     *
     * Performs validation on the band structure and pairing.
     * Combines both left and right scans to find all base pairs and their children.
     *
     *
     * @param lb Left border position.
     * @param li Left inner position.
     * @param ri Right inner position.
     * @param rb Right border position.
     * @param pair_table Base-pair index mapping for the structure.
     * @param node_table Vector of pointers to LoopNode objects for each position in the structure.
     * @throws DetailedException if band structure is invalid.
     */
    static Band construct_band(std::size_t lb, std::size_t li, std::size_t ri, std::size_t rb,
                               const std::vector<std::size_t>& pair_table,
                               const std::vector<LoopNode*>& node_table);

    // Convenience method (BandBounds stores lb, li, ri, rb)
    static Band construct_band(BandBounds bounds, const std::vector<std::size_t>& pair_table,
                               const std::vector<LoopNode*>& node_table);

   private:
    /**
     * @brief Build the list of base pairs that participate in the band.
     *
     * Scans left arm for base pairs that are part of the band.
     * Stores nested children only from the left arm.
     *
     * @param lb Left border position.
     * @param li Left inner position.
     * @param ri Right inner position.
     * @param rb Right border position.
     * @param pair_table Base-pair index mapping for the structure.
     * @param child_count Reference to the count of child nodes.
     * @return Vector of PKBasePair objects representing the base pairs in the band.
     *
     */
    static std::vector<PKBasePair> find_base_pairs_left_scan(
        std::size_t lb, std::size_t li, std::size_t ri, std::size_t rb,
        const std::vector<std::size_t>& pair_table, const std::vector<LoopNode*>& node_table,
        int& child_count);

    /**
     * @brief Populate closed region children from children of the right arm.
     *
     * @param base_pairs Vector of PKBasePair objects to populate.
     * @param ri Right inner position.
     * @param rb Right border position.
     * @param node_table Vector of pointers to LoopNode objects for each position in the structure.
     * @param child_count Reference to the count of child nodes.
     */
    static void populate_right_arm_children(std::vector<PKBasePair>& base_pairs, std::size_t ri,
                                            std::size_t rb,
                                            const std::vector<LoopNode*>& node_table,
                                            int& child_count);
};

}  // namespace knotergy