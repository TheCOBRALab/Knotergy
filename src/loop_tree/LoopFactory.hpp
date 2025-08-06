#pragma once

#include "../preprocessing/RNAProcessedEntry.hpp"
#include "LoopNode.hpp"

namespace knotergy {
class LoopFactory {
   public:
    LoopFactory(const RNAProcessedEntry& processed_rna);

    std::shared_ptr<LoopNode> get_root_node() { return root_node_; };

    void print_tree(bool debug = false) const;
    void print_tree(const std::shared_ptr<LoopNode>& node, size_t depth, bool debug = false) const;

   private:
    const RNAProcessedEntry& processed_rna_;
    std::shared_ptr<LoopNode> root_node_;
    size_t structure_length_;

    /**
     * @brief Constructs a hierarchical tree of loop regions from a list of closed regions.
     *
     * This method builds a loop tree where each node represents a closed region (a contiguous
     * segment of a secondary structure enclosed by base pairs). It uses a stack-based approach to
     * assign parent-child relationships based on nesting (i.e., a region fully contained in another
     * becomes its child).
     *
     * The root node represents the external loop (entire structure), and all other regions are
     * nested inside it according to their boundaries. Loop types, unpaired base counts, band
     * annotations, and pseudo-nested structures are also determined in the process.
     *
     * @param closed_regions A list of closed regions to be organized into a loop tree.
     *
     * @details
     * - Regions are first bucket-sorted by their starting indices.
     * - A root node is created to span the entire sequence range.
     * - A stack is used to maintain the current active path in the tree as nodes are added.
     * - Each region is added in order, and completed nodes are processed for loop type,
     *   band annotations, and pseudo-nesting checks.
     * - The method finalizes remaining nodes on the stack after all regions have been added.
     *
     * @note
     * - Assumes `processed_rna` provides access to the secondary structure, base pairings, and
     * unpaired counts.
     * - `NULL_INDEX` is used to represent an invalid or out-of-bound index for the root node.
     * - The tree is stored in `root_node_`, and child nodes are recursively linked.
     */
    void build_tree(const std::vector<ClosedRegion>& closed_regions);

    void populate_node(LoopNode& node);
    void populate_node(const std::shared_ptr<LoopNode>& node);

    int count_unpaired_bases_excluding_children(const LoopNode& node);
    LoopType find_loop_type(const LoopNode& node);
    void pseudo_nested_check(LoopNode& node);
    void label_pseudonested_children(LoopNode& node);
    std::vector<ClosedRegion> closed_region_bucket_sort(
        const std::vector<ClosedRegion>& closed_regions);
    void annotate_bands(const std::shared_ptr<LoopNode>& node);
};

}  // namespace knotergy