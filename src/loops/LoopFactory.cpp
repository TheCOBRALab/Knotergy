#include "LoopFactory.hpp"

#include <algorithm>
#include <memory>
#include <stack>

#include "../loops/BandFinder.hpp"

namespace knotergy {

LoopFactory::LoopFactory(const RNAProcessedEntry& processed_rna) : processed_rna_(processed_rna) {
    build_tree(processed_rna.get_closed_regions());
}

// Builds a tree of closed regions. Each node represents a closed region
// Children of a node are all stems directly nested inside the closed region
void LoopFactory::build_tree(const std::vector<ClosedRegion>& closed_regions) {
    // Sort regions by starting index (ascending)
    size_t structure_length = processed_rna_.size();
    std::vector<ClosedRegion> sorted_closed_regions =
        closed_region_bucket_sort(closed_regions, structure_length);

    // Root node covers full structure [-1, structure_length]
    // Use NULL_INDEX for -1 (unsigned type)
    // Indicies are out of bounds so the root node is always larger than any closed region inside
    root_node_ = std::make_shared<LoopNode>(ClosedRegion{NULL_INDEX, structure_length});
    root_node_->loop_type = LoopType::External;

    std::stack<std::shared_ptr<LoopNode>> node_stack;
    node_stack.push(root_node_);

    BandFinder band_finder(processed_rna_.get_pairings());

    // Pop until node_stack.end() is the parent of current node
    // A node is only popped (and processed) after all of its children have been added.
    // Therefore, its loop type, bands and pseudo-nested bands can now be determined.
    for (const ClosedRegion& closed_region : sorted_closed_regions) {
        while (node_stack.top()->end < closed_region.begin) {
            std::shared_ptr<LoopNode>& node = node_stack.top();
            node->loop_type = find_loop_type(*node);
            band_finder.annotate_bands(node);
            count_unpaired_bases_excluding_children(*node);
            pseudo_nested_check(*node);
            node_stack.pop();
        }

        // Parent = parent of current node. Child = current node
        // Creates child node and links it to its parent
        std::shared_ptr<LoopNode>& parent = node_stack.top();
        parent->children.emplace_back(std::make_shared<LoopNode>(closed_region));
        std::shared_ptr<LoopNode>& child = parent->children.back();
        child->parent = parent;  // May be unused

        // Gets the total unpaired base pairs within the closed region (including that of children)
        child->total_number_of_unpaired_bases = processed_rna_.get_unpaired_count(closed_region);

        node_stack.push(child);
    }

    // process all remaining nodes
    while (node_stack.top()->begin != NULL_INDEX) {
        std::shared_ptr<LoopNode>& node = node_stack.top();
        count_unpaired_bases_excluding_children(*node);
        node->loop_type = find_loop_type(*node);
        band_finder.annotate_bands(node);
        pseudo_nested_check(*node);
        node_stack.pop();
    }
}

void LoopFactory::count_unpaired_bases_excluding_children(LoopNode& node) {
    int total = node.total_number_of_unpaired_bases;
    for (std::shared_ptr<LoopNode> child : node.children) {
        total -= child->total_number_of_unpaired_bases;
    }
    node.number_of_exclusive_unpaired_bases = total;
}

LoopType LoopFactory::find_loop_type(const LoopNode& node) {
    const std::vector<size_t>& pairings = processed_rna_.get_pairings();

    if (pairings[node.begin] != node.end) {
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
        // ...and if not stacked, then bulge or internal loop
        return LoopType::Internal;
    }

    // Multiple children or a pseudoknotted child = Multiloop
    return LoopType::Multibranch;
}

void LoopFactory::pseudo_nested_check(LoopNode& node) {
    if (node.loop_type != LoopType::Pseudoknot) return;

    for (std::shared_ptr<LoopNode> child_node : node.children) {
        if (child_node->pseudo_type == PseudoNestedType::InsideBand) {
            ++node.number_of_children_inside_band;
        } else if (child_node->pseudo_type == PseudoNestedType::OutsideBand) {
            ++node.number_of_children_outside_band;
            node.number_of_unpaired_bases_in_children_outside_band +=
                static_cast<int>(child_node->end - child_node->begin + 1);
        }
    }
}

std::vector<ClosedRegion> LoopFactory::closed_region_bucket_sort(
    const std::vector<ClosedRegion>& closed_regions, size_t structure_length) {
    /*  one linear pass buckets the regions by their left index  */
    std::vector<std::vector<ClosedRegion>> buckets(structure_length);
    for (const ClosedRegion& cr : closed_regions) {
        if (!buckets[cr.begin].empty()) {
            throw std::runtime_error("Duplicate starting index `" + std::to_string(cr.begin) +
                                     "` found in closed regions");
        }
        buckets[cr.begin].push_back(cr);
    }

    /*  then walk the buckets left-to-right: O(n)  */
    std::vector<ClosedRegion> sorted_regions;
    sorted_regions.reserve(closed_regions.size());
    for (std::vector<ClosedRegion> b : buckets) {
        for (ClosedRegion& cr : b) {
            sorted_regions.push_back(cr);
        }
    }

    return sorted_regions;
}

void LoopFactory::print_tree(bool debug) const {
    for (const auto& child : root_node_->children) {
        print_tree(child, 0, debug);
    }
}

void LoopFactory::print_tree(const std::shared_ptr<LoopNode>& node, size_t depth,
                             bool debug) const {
    std::cout << std::string(depth, '.')  // indent with dots
              << '[' << node->begin << ',' << node->end << "]  "
              << "  unpaired=" << node->number_of_exclusive_unpaired_bases
              << "  children=" << node->children.size() << '\n';

    if (debug) {
        std::cout << *node << '\n';
    }

    for (const auto& child : node->children) {
        print_tree(child, depth + 1, debug);
    }
}

}  // namespace knotergy