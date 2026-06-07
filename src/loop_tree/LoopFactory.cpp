#include "LoopFactory.hpp"

#include <algorithm>
#include <memory>
#include <stack>
#include <unordered_set>

#include "../loop_tree/BandFinder.hpp"

namespace knotergy {

LoopFactory::LoopFactory(const ProcessedRNAEntry& processed_rna) : processed_rna_{processed_rna} {
    // CLOSED REGIONS MUST BE SORTED BY START INDEX FOR THE BUILDING ALGORITHM TO WORK CORRECTLY
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

    root_node_ = std::make_unique<LoopNode>(ClosedRegion{NULL_INDEX, processed_rna_.size()});
    root_node_->loop_type = LoopType::External;

    std::vector<LoopNode*> node_stack;
    node_stack.reserve(closed_regions.size() + 1);
    node_stack.push_back(root_node_.get());

    for (const ClosedRegion& cr : closed_regions) {
        // Pop until node_stack.end() is the parent of current node
        // A node is only popped (and processed) after all of its children have been added.
        // Therefore, its loop type, bands and pseudo-nested bands can now be determined.
        while (node_stack.back()->end < cr.begin) {
            populate_node(*node_stack.back());
            node_stack.pop_back();
        }

        // parent = parent of current node. child = current node
        LoopNode* parent = node_stack.back();

        std::unique_ptr<LoopNode> child = std::make_unique<LoopNode>(cr);
        child->parent = parent;
        child->total_unpaired_bases_count = processed_rna_.get_unpaired_count(cr);

        LoopNode* child_raw = child.get();
        parent->children.emplace_back(std::move(child));

        node_stack.push_back(child_raw);
    }

    // process all remaining nodes
    while (node_stack.back()->begin != NULL_INDEX) {
        LoopNode* node = node_stack.back();
        populate_node(*node);
        node_stack.pop_back();
    }
}

void LoopFactory::populate_node(LoopNode& node) {
    node.exclusive_unpaired_bases_count = count_unpaired_bases_excluding_children(node);
    node.loop_type = find_loop_type(node);

    if (node.loop_type == LoopType::Pseudoknot) {
        // Used for band finding and navigation.
        // Lazily initialize when we encounter the first pseudoknot.
        if (aux_bands_.empty()) {
            aux_bands_.resize(processed_rna_.get_structure().size());
        }

        node.bands = BandFinder::find_bands(node, aux_bands_, processed_rna_);
        node.total_number_of_base_pairs = count_total_base_pairs(node);
        label_pseudonested_children(node);
        pseudo_nested_check(node);
    }
}

void LoopFactory::populate_node(const std::unique_ptr<LoopNode>& node) {
    populate_node(*node);
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
    for (const std::unique_ptr<LoopNode>& child : node.children) {
        total -= child->total_unpaired_bases_count;
    }
    return total;
}

LoopType LoopFactory::find_loop_type(const LoopNode& node) {
    const std::vector<size_t>& pair_table = processed_rna_.get_pair_table();

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
        // ...and if not stacked, then bulge or internal loop
        return LoopType::Internal;
    }

    // Multiple children or a pseudoknotted child = Multiloop
    return LoopType::Multibranch;
}

void LoopFactory::pseudo_nested_check(LoopNode& node) {
    if (node.loop_type != LoopType::Pseudoknot) return;

    for (const std::unique_ptr<LoopNode>& child_node : node.children) {
        if (child_node->pseudo_type == PseudoNestedType::WithinBand) {
            ++node.number_of_withinband_children;
        } else if (child_node->pseudo_type == PseudoNestedType::Nested) {
            ++node.number_of_nested_children;
        }
    }
}

// Expects bands to already be populated for the node
void LoopFactory::label_pseudonested_children(LoopNode& node) {
    /* only pseudoknots need bands; leave the rest untouched */
    if (node.loop_type != LoopType::Pseudoknot) return;
    if (node.children.empty()) return;

    /**
     * There are two methods: Linear & Non-Linear
     *
     * Non-Linear is usually more efficient, but in cases where there are many bands and
     * many children, the linear method can be more efficient.
     */

    // Non-linear method is preferred when it has less operations than the linear method.
    if (static_cast<size_t>(node.total_number_of_base_pairs) >=
        node.children.size() * node.bands.size()) {
        for (const std::unique_ptr<LoopNode>& child_node : node.children) {
            child_node->pseudo_type = PseudoNestedType::Nested;
        }
        for (const Band& band : node.bands) {
            if (band.get_number_of_children() == 0) continue;
            for (const std::unique_ptr<LoopNode>& child_node : node.children) {
                // If start index is in the band, by definition, the end index must also be in the
                // band.
                if (child_node && band.contains(child_node->begin)) {
                    child_node->pseudo_type = PseudoNestedType::WithinBand;
                };
            }
        }
    } else {  // Linear method is preferred when it has less operations than the non-linear method.
        std::unordered_set<size_t> within_band_start_idx;
        within_band_start_idx.reserve(node.children.size());

        for (const Band& band : node.bands) {
            if (band.get_number_of_children() == 0) continue;
            for (const BasePair& base_pair : band.base_pairs()) {
                within_band_start_idx.insert(base_pair.i);
            }
        }
        for (const std::unique_ptr<LoopNode>& child_node : node.children) {
            if (within_band_start_idx.count(child_node->begin)) {
                child_node->pseudo_type = PseudoNestedType::WithinBand;
            } else {
                child_node->pseudo_type = PseudoNestedType::Nested;
            }
        }
    }
}

void LoopFactory::print_tree(bool debug) const {
    for (const auto& child : root_node_->children) {
        print_tree(child, 0, debug);
    }
}

void LoopFactory::print_tree(const std::unique_ptr<LoopNode>& node, size_t depth,
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

    std::vector<std::unique_ptr<LoopNode>> work;
    work.emplace_back(std::move(root_node_));

    while (!work.empty()) {
        auto node = std::move(work.back());
        work.pop_back();

        for (auto& child : node->children) {
            if (child) work.emplace_back(std::move(child));
        }

        node->children.clear();
        // no parent reset needed; parent is now a raw non-owning pointer
    }
}

}  // namespace knotergy