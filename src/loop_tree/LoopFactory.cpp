#include "LoopFactory.hpp"

#include <algorithm>
#include <memory>
#include <stack>

namespace knotergy {

LoopFactory::LoopFactory(const ProcessedRNAEntry& processed_rna) : pRNA_{processed_rna} {
    node_table_.resize(pRNA_.size(), nullptr);

    // Reserve space for all nodes (including root node)
    // NOTICE: Avoid reallocations. That would invalidate pointers.
    nodes_.reserve(pRNA_.get_closed_regions().size() + 1);  // +1 for root node

    // CLOSED REGIONS MUST BE SORTED BY START INDEX FOR THE BUILDING ALGORITHM TO WORK CORRECTLY
    // A check is performed in build_tree() to ensure this precondition is met.
    // pre populate node_table_ with nullptrs to reserve space for all nodes
    build_tree(processed_rna.get_closed_regions());
}

// Builds a tree of closed regions. Each node represents a closed region
// Children of a node are all stems directly nested inside the closed region
void LoopFactory::build_tree(const std::vector<ClosedRegion>& closed_regions) {
    // Root node covers full structure [-1, structure_length]
    // Use NULL_INDEX for -1 (unsigned type)

    // This should never happen because RNAProcessor::compute_closed_regions() guarantees sorted
    // output, but we check just in case to avoid silently producing incorrect loop trees.
    if (!std::is_sorted(closed_regions.begin(), closed_regions.end())) {
        THROW_ERROR(
            "Closed regions are not sorted by start index. This may cause the loop factory "
            "algorithm to fail.\n");
    }

    nodes_.emplace_back(ClosedRegion{NULL_INDEX, pRNA_.size()});
    root_node_ = &nodes_.back();
    root_node_->loop_type = LoopType::External;

    std::vector<LoopNode*> node_stack;
    node_stack.reserve(closed_regions.size() + 1);
    node_stack.push_back(root_node_);

    for (const ClosedRegion& cr : closed_regions) {
        // Pop until node_stack.end() is the parent of current node
        // A node is only popped (and processed) after all of its children have been added.
        // Therefore, its loop type, bands and pseudo-nested bands can now be determined.
        while (node_stack.back()->end < cr.begin) {
            LoopNode* node = node_stack.back();

            node_table_[node->begin] = node;
            node_table_[node->end] = node;

            populate_node(*node);
            node_stack.pop_back();
        }

        // parent = parent of current node. child = current node
        LoopNode* parent = node_stack.back();
        nodes_.emplace_back(cr);

        LoopNode* child = &nodes_.back();
        child->parent = parent;
        child->total_unpaired_bases_count = pRNA_.get_unpaired_count(cr);

        parent->children.emplace_back(child);

        node_stack.push_back(child);
    }

    // process all remaining nodes
    while (node_stack.back()->begin != NULL_INDEX) {
        LoopNode* node = node_stack.back();

        node_table_[node->begin] = node;
        node_table_[node->end] = node;

        populate_node(*node);
        node_stack.pop_back();
    }
}

void LoopFactory::populate_node(LoopNode& node) {
    node.exclusive_unpaired_bases_count = count_unpaired_bases_excluding_children(node);
    node.loop_type = find_loop_type(node);
    // Populate node encodings for ViennaRNA energy calculations.

    if (node.loop_type == LoopType::Pseudoknot) {
        // Used for band finding and navigation.
        // Lazily initialize when we encounter the first pseudoknot.
        if (aux_bands_.empty()) {
            aux_bands_.resize(pRNA_.get_structure().size());
        }

        node.bands = BandFinder::find_bands(node, aux_bands_, pRNA_, node_table_);
        node.total_number_of_base_pairs = count_total_base_pairs(node);
        pseudo_nested_check(node);
    }
}

int LoopFactory::count_total_base_pairs(const LoopNode& node) {
    size_t total = 0;
    for (const Band& band : node.bands) {
        total += band.base_pairs().size();
    }
    return static_cast<int>(total);
}

int LoopFactory::count_unpaired_bases_excluding_children(const LoopNode& node) {
    int total = node.total_unpaired_bases_count;
    for (const LoopNode* child : node.children) {
        total -= child->total_unpaired_bases_count;
    }
    return total;
}

LoopType LoopFactory::find_loop_type(const LoopNode& node) {
    const std::vector<size_t>& pair_table = pRNA_.get_pair_table();

    if (pair_table[node.begin] != node.end) {
        return LoopType::Pseudoknot;
    }

    // Zero Children
    if (node.children.empty()) {
        return LoopType::Hairpin;
    }

    // if only one child and it's not pseudoknotted...
    if (node.children.size() == 1 && node.children[0]->loop_type != LoopType::Pseudoknot) {
        // ... and if child's base pair is stacked with this loop's pair, then stacked
        if ((node.children[0]->begin == node.begin + 1) && node.children[0]->end == node.end - 1) {
            return LoopType::Stack;
        }
        // ...and if not stacked, and the child is adjacent to either border, then it's a bulge
        if (node.children[0]->begin == node.begin + 1 || node.children[0]->end == node.end - 1) {
            return LoopType::Bulge;
        }
        // ...otherwise, it's an internal loop.
        return LoopType::Internal;
    }

    // Multiple children or a pseudoknotted child = Multiloop
    return LoopType::Multibranch;
}

void LoopFactory::pseudo_nested_check(LoopNode& node) {
    if (node.loop_type != LoopType::Pseudoknot) return;
    if (node.children.empty()) return;

    // Count the number of children that are within bands and outside bands.
    if (node.bands.size() <= node.children.size()) {
        for (const Band& band : node.bands) {
            node.number_of_withinband_children += band.get_number_of_children();
        }
        int total_children = static_cast<int>(node.children.size());
        node.number_of_outsideband_children = total_children - node.number_of_withinband_children;
    } else {
        for (const LoopNode* child_node : node.children) {
            if (child_node->pseudo_type == PseudoNestedType::WithinBand) {
                ++node.number_of_withinband_children;
            } else if (child_node->pseudo_type == PseudoNestedType::OutsideBandIntervals) {
                ++node.number_of_outsideband_children;
            }
        }
    }
}

void LoopFactory::print_tree(bool debug) const {
    for (const LoopNode* child : root_node_->children) {
        print_tree(child, 0, debug);
    }
}

void LoopFactory::print_tree(const LoopNode* node, size_t depth, bool debug) const {
    std::cout << std::string(depth, '.')  // indent with dots
              << '[' << node->begin << ',' << node->end << "]  "
              << "  unpaired=" << node->exclusive_unpaired_bases_count
              << "  children=" << node->children.size() << '\n';

    if (debug) {
        std::cout << *node << '\n';
    }

    for (const LoopNode* child : node->children) {
        print_tree(child, depth + 1, debug);
    }
}

}  // namespace knotergy