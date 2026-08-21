#pragma once

#include "energy/vienna/ViennaUtils.hpp"
#include "loop_tree/LoopNode.hpp"
#include "loop_tree/bands/BandFinder.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

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
    ~LoopFactory() { destroy_tree_iterative(); };

    /**
     * @brief Get the root node of the loop tree.
     *
     * The root node represents the external loop containing all other loops.
     *
     * @return Raw pointer to the root LoopNode.
     */
    [[nodiscard]] LoopNode& get_root_node() { return *root_node_; }

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
    void print_tree(const std::unique_ptr<LoopNode>& node, size_t depth, bool debug = false) const;

   private:
    const ProcessedRNAEntry& pRNA_;
    std::unique_ptr<LoopNode> root_node_;
    size_t structure_length_;
    std::vector<PairedBaseNode> aux_bands_;  ///< Auxiliary structure for band finder.
    std::vector<LoopNode*>
        node_table_;  ///< All nodes for closed region skipping and pseudo-nested checks.

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
     * - Assumes `processed_rna` provides access to the secondary structure, base pair_table, and
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
     * @brief Populate a loop node (unique pointer version).
     *
     * @param node Unique pointer to the loop node to populate.
     */
    void populate_node(const std::unique_ptr<LoopNode>& node);

    /**
     * @brief Count the total number of base pairs in a loop node.
     *
     * @param node The loop node to analyze.
     * @return Total number of base pairs in the loop's closed region (excluding children).
     */
    [[nodiscard]] int count_total_base_pairs(const LoopNode& node);

    /**
     * @brief Count unpaired bases in a loop, excluding those in child loops.
     *
     * @param node The loop node to analyze.
     * @return Number of unpaired bases exclusive to this loop.
     */
    [[nodiscard]] int count_unpaired_bases_excluding_children(const LoopNode& node);

    /**
     * @brief Determine the loop type of a node.
     *
     * Classifies the loop as Stack, Hairpin, Internal, Multibranch, External, or Pseudoknot.
     *
     * @param node The loop node to classify.
     * @return The determined LoopType.
     */
    [[nodiscard]] LoopType find_loop_type(const LoopNode& node);

    /**
     * @brief Check if a loop has pseudo-nested children and update pseudo_type.
     *
     * @param node The loop node to check (modified in place).
     */
    void pseudo_nested_check(LoopNode& node);

    /**
     * @brief Annotate a loop node with pseudoknot band information.
     *
     * Identifies and stores Band objects for pseudoknotted loops.
     *
     * @param node Unique pointer to the loop node to annotate.
     */
    void annotate_bands(const std::unique_ptr<LoopNode>& node);

    /**
     * @brief Iteratively destroy the loop tree to free memory.
     *
     * This method uses a stack to traverse the tree and delete nodes without recursion,
     * preventing potential stack overflow for deeply nested structures.
     */
    void destroy_tree_iterative();
};

}  // namespace knotergy