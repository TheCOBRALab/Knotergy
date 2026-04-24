#pragma once

#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "LoopNode.hpp"
#include "BandFinder.hpp"

namespace knotergy {
/**
 * @brief Factory class for building loop tree representations of RNA secondary structures.
 *
 * LoopFactory constructs a hierarchical tree of LoopNode objects from a ProcessedRNAEntry.
 * Each node represents a loop region (external, hairpin, internal, multibranch, or pseudoknot).
 * The tree structure captures the nesting relationships between loops.
 */
class LoopFactory {
   public:
    /**
     * @brief Construct a LoopFactory and build the loop tree.
     *
     * Upon construction, immediately builds the complete loop tree from the
     * processed RNA entry's closed regions.
     *
     * @param processed_rna The processed RNA entry containing structure and pairing information.
     */
    LoopFactory(const ProcessedRNAEntry& processed_rna);
    
    /**
     * @brief Destroy the LoopFactory and its associated loop tree.
     * 
     * This destructor ensures that all dynamically allocated LoopNode objects are properly
     * deallocated to prevent memory leaks. This is needed for deeply nested structures.
     */
    ~LoopFactory() {destroy_tree_iterative();};

    /**
     * @brief Get the root node of the loop tree.
     *
     * The root node represents the external loop containing all other loops.
     *
     * @return Shared pointer to the root LoopNode.
     */
    std::shared_ptr<LoopNode> get_root_node() { return root_node_; };

    /**
     * @brief Print the loop tree structure.
     *
     * @param debug Whether to include debug information in the output (default: false).
     */
    void print_tree(bool debug = false) const;

    /**
     * @brief Print a loop tree node and its children recursively.
     *
     * @param node The node to print.
     * @param depth Current depth in the tree (for indentation).
     * @param debug Whether to include debug information (default: false).
     */
    void print_tree(const std::shared_ptr<LoopNode>& node, size_t depth, bool debug = false) const;

    /**
     * @brief Bucket sort closed regions by their starting position.
     *
     * @param closed_regions Vector of closed regions to sort.
     * @param structure_length Length of the RNA structure.
     * @return Sorted vector of closed regions.
     */
    static std::vector<ClosedRegion> closed_region_bucket_sort(
        const std::vector<ClosedRegion>& closed_regions, size_t structure_length);

   private:
    const ProcessedRNAEntry& processed_rna_;
    std::shared_ptr<LoopNode> root_node_;
    size_t structure_length_;
    std::vector<PairedBaseNode> aux_bands_;  ///< Auxiliary structure for band detection and navigation.

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

    /**
     * @brief Populate a loop node with type, unpaired counts, and band information.
     *
     * @param node The loop node to populate.
     */
    void populate_node(LoopNode& node);

    /**
     * @brief Populate a loop node (shared pointer version).
     *
     * @param node Shared pointer to the loop node to populate.
     */
    void populate_node(const std::shared_ptr<LoopNode>& node);

    /**
     * @brief Count unpaired bases in a loop, excluding those in child loops.
     *
     * @param node The loop node to analyze.
     * @return Number of unpaired bases exclusive to this loop.
     */
    int count_unpaired_bases_excluding_children(const LoopNode& node);

    /**
     * @brief Determine the loop type of a node.
     *
     * Classifies the loop as Stack, Hairpin, Internal, Multibranch, External, or Pseudoknot.
     *
     * @param node The loop node to classify.
     * @return The determined LoopType.
     */
    LoopType find_loop_type(const LoopNode& node);

    /**
     * @brief Check if a loop has pseudo-nested children and update pseudo_type.
     *
     * @param node The loop node to check (modified in place).
     */
    void pseudo_nested_check(LoopNode& node);

    /**
     * @brief Label child nodes that are pseudo-nested or within bands.
     *
     * @param node The parent loop node containing children to label.
     */
    void label_pseudonested_children(LoopNode& node);

    /**
     * @brief Annotate a loop node with pseudoknot band information.
     *
     * Identifies and stores Band objects for pseudoknotted loops.
     *
     * @param node Shared pointer to the loop node to annotate.
     */
    void annotate_bands(const std::shared_ptr<LoopNode>& node);

    /**
     * @brief Iteratively destroy the loop tree to free memory.
     *
     * This method uses a stack to traverse the tree and delete nodes without recursion,
     * preventing potential stack overflow for deeply nested structures.
     */
    void destroy_tree_iterative();
};

}  // namespace knotergy