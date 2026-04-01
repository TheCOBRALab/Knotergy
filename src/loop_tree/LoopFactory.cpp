#include "LoopFactory.hpp"

#include <algorithm>
#include <memory>
#include <stack>
#include <unordered_set>

#include "../loop_tree/BandFinder.hpp"

namespace knotergy {

LoopFactory::LoopFactory(const ProcessedRNAEntry& processed_rna) : processed_rna_{processed_rna} {
    build_tree(processed_rna.get_closed_regions());
}

// Builds a tree of closed regions. Each node represents a closed region
// Children of a node are all stems directly nested inside the closed region
void LoopFactory::build_tree(const std::vector<ClosedRegion>& closed_regions) {
    // Sort regions by starting index (ascending)
    std::vector<ClosedRegion> sorted_closed_regions =
        closed_region_bucket_sort(closed_regions, processed_rna_.size());

    // Root node covers full structure [-1, structure_length]
    // Use NULL_INDEX for -1 (unsigned type)
    root_node_ = std::make_shared<LoopNode>(ClosedRegion{NULL_INDEX, processed_rna_.size()});
    root_node_->loop_type = LoopType::External;

    std::stack<std::shared_ptr<LoopNode>> node_stack;
    node_stack.push(root_node_);

    for (const ClosedRegion& closed_region : sorted_closed_regions) {
        // Pop until node_stack.end() is the parent of current node
        // A node is only popped (and processed) after all of its children have been added.
        // Therefore, its loop type, bands and pseudo-nested bands can now be determined.
        while (node_stack.top()->end < closed_region.begin) {
            std::shared_ptr<LoopNode> node = node_stack.top();
            populate_node(node);
            node_stack.pop();
        }

        // parent = parent of current node. child = current node
        std::shared_ptr<LoopNode> parent = node_stack.top();
        std::shared_ptr<LoopNode> child = std::make_shared<LoopNode>(closed_region);
        child->parent = parent;
        parent->children.emplace_back(child);

        // Gets the total unpaired base pairs within the closed region (including that of children)
        child->total_unpaired_bases_count = processed_rna_.get_unpaired_count(closed_region);

        node_stack.push(child);
    }

    // process all remaining nodes
    while (node_stack.top()->begin != NULL_INDEX) {
        std::shared_ptr<LoopNode>& node = node_stack.top();
        populate_node(node);
        node_stack.pop();
    }
}

void LoopFactory::populate_node(LoopNode& node) {
    node.exclusive_unpaired_bases_count = count_unpaired_bases_excluding_children(node);
    node.loop_type = find_loop_type(node);
    node.bands = BandFinder::find_bands(node, processed_rna_);
    node.number_of_bands = static_cast<int>(node.bands.size());
    label_pseudonested_children(node);
    pseudo_nested_check(node);
}

void LoopFactory::populate_node(const std::shared_ptr<LoopNode>& node) {
    populate_node(*node);
}

int LoopFactory::count_unpaired_bases_excluding_children(const LoopNode& node) {
    int total = node.total_unpaired_bases_count;
    for (std::shared_ptr<LoopNode> child : node.children) {
        total -= child->total_unpaired_bases_count;
    }
    return total;
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
        if (child_node->pseudo_type == PseudoNestedType::WithinBand) {
            ++node.number_of_withinband_children;
        } else if (child_node->pseudo_type == PseudoNestedType::Nested) {
            ++node.number_of_nested_children;
        }
    }
}

void LoopFactory::label_pseudonested_children(LoopNode& node) {
    /* only pseudoknots need bands; leave the rest untouched */
    if (node.loop_type != LoopType::Pseudoknot) return;

    const std::vector<size_t>& cr_pairings = processed_rna_.get_closed_regions_pairings();
    std::unordered_set<size_t> within_band_start_idx;
    within_band_start_idx.reserve(node.children.size());

    for (const Band& band : node.bands) {
        for (size_t i = band.left_border() + 1; i <= band.left_inner(); ++i) {
            if (cr_pairings[i] != NULL_INDEX && (i < cr_pairings[i])) {
                within_band_start_idx.emplace(i);
                i = cr_pairings[i];
                continue;
            }
        }
        for (size_t i = band.right_inner() + 1; i <= band.right_border(); ++i) {
            if (cr_pairings[i] != NULL_INDEX && (i < cr_pairings[i])) {
                within_band_start_idx.emplace(i);
                i = cr_pairings[i];
                continue;
            }
        }
    }

    for (std::shared_ptr<LoopNode> c : node.children) {
        if (within_band_start_idx.find(c->begin) != within_band_start_idx.end()) {
            c->pseudo_type = PseudoNestedType::WithinBand;
        } else {
            c->pseudo_type = PseudoNestedType::Nested;
        }
    }
}

std::vector<ClosedRegion> LoopFactory::closed_region_bucket_sort(
    const std::vector<ClosedRegion>& closed_regions, const size_t& structure_length) {
    // Map each possible 'begin' position to an index in closed_regions.
    std::vector<size_t> begin_to_idx(structure_length, NULL_INDEX);

    for (size_t i = 0; i < closed_regions.size(); ++i) {
        size_t begin = closed_regions[i].begin;
        if (begin >= structure_length) {
            THROW_ERROR("Begin `" + std::to_string(begin) + "` out of range [0," +
                        std::to_string(structure_length) + ")");
        }
        if (begin_to_idx[begin] != SIZE_MAX) {
            THROW_ERROR("Duplicate starting index `" + std::to_string(begin) +
                        "` found in closed regions");
        }
        begin_to_idx[begin] = i;
    }

    std::vector<ClosedRegion> sorted;
    sorted.reserve(closed_regions.size());
    for (size_t begin = 0; begin < structure_length; ++begin) {
        const size_t idx = begin_to_idx[begin];
        if (idx != SIZE_MAX) sorted.push_back(closed_regions[idx]);
    }
    return sorted;
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
              << "  unpaired=" << node->exclusive_unpaired_bases_count
              << "  children=" << node->children.size() << '\n';

    if (debug) {
        std::cout << *node << '\n';
    }

    for (const auto& child : node->children) {
        print_tree(child, depth + 1, debug);
    }
}

// Iteratively destroy the loop tree to free memory.
void LoopFactory::destroy_tree_iterative() {
        if (!root_node_) return;

        std::vector<std::shared_ptr<LoopNode>> work;
        work.push_back(std::move(root_node_));

        while (!work.empty()) {
            auto node = std::move(work.back());
            work.pop_back();

            for (auto& child : node->children) {
                if (child) work.push_back(std::move(child));
            }

            node->children.clear();
            node->parent.reset();
        }
    }

}  // namespace knotergy