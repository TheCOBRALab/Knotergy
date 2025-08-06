#include "LoopFactory.hpp"

#include <algorithm>
#include <memory>
#include <stack>

#include "../loop_tree/BandFinder.hpp"

namespace knotergy {

LoopFactory::LoopFactory(const RNAProcessedEntry& processed_rna) : processed_rna_{processed_rna} {
    build_tree(processed_rna.get_closed_regions());
}

// Builds a tree of closed regions. Each node represents a closed region
// Children of a node are all stems directly nested inside the closed region
void LoopFactory::build_tree(const std::vector<ClosedRegion>& closed_regions) {
    // Sort regions by starting index (ascending)
    std::vector<ClosedRegion> sorted_closed_regions = closed_region_bucket_sort(closed_regions);

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
            std::shared_ptr<LoopNode>& node = node_stack.top();
            populate_node(node);
            node_stack.pop();
        }

        // parent = parent of current node. child = current node
        std::shared_ptr<LoopNode>& parent = node_stack.top();
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
    node.bands = BandFinder::find_bands(node, processed_rna_.get_pairings(), processed_rna_.get_closed_regions_pairings());
    node.number_of_bands = static_cast<int>(node.bands.size());
    label_pseudonested_children(node);
    pseudo_nested_check(node);
}

void LoopFactory::populate_node(const std::shared_ptr<LoopNode>& node){
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
        if (child_node->pseudo_type == PseudoNestedType::InsideBand) {
            ++node.number_of_insideband_children;
        } else if (child_node->pseudo_type == PseudoNestedType::CrossBand) {
            ++node.number_of_crossband_children;
            node.number_of_unpaired_bases_in_crossband_children +=
                static_cast<int>(child_node->end - child_node->begin - 1);
        }
    }
}

void LoopFactory::label_pseudonested_children(LoopNode& node) {
        /* only pseudoknots need bands; leave the rest untouched */
        if (node.loop_type != LoopType::Pseudoknot) return;

        const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
        const std::vector<Band>& bands = node.bands;
        size_t child_idx = 0;
        size_t band_idx = 0;
        
        // checks if a child is exclusively in one band (InsideBand), or if its CrossBand 
        while (child_idx < children.size()) {
            std::shared_ptr<LoopNode> child = children[child_idx];
            child->pseudo_type = PseudoNestedType::CrossBand;

            // Prevent out-of-bounds in band lookup
            if (band_idx >= bands.size()) {
                ++child_idx;
                continue;
            }

            // Checks if child is exclusively in one band
            if (node.bands[band_idx].contains(child->begin, child->end)) {
                bool crosses_previous = 
                (band_idx > 0) && bands[band_idx - 1].contains(child->begin, child->end);
 
                bool crosses_next = 
                (band_idx + 1) < bands.size() && bands[band_idx + 1].contains(child->begin, child->end);

                if (!crosses_previous && !crosses_next) {
                    child->pseudo_type = PseudoNestedType::InsideBand;
                }
                ++child_idx;  // Done with this child, move to next
            } else {
                ++band_idx;
            }
        }
    }

std::vector<ClosedRegion> LoopFactory::closed_region_bucket_sort(const std::vector<ClosedRegion>& closed_regions) {
    // Map regions by their starting index
    // Reserve capacity to avoid rehashing
    std::unordered_map<size_t, ClosedRegion> index_to_region;
    index_to_region.reserve(closed_regions.size());

    size_t min_index = NULL_INDEX;
    size_t max_index = 0;

    for (const ClosedRegion& cr : closed_regions) {
        // Try inserting; .second is false if this start index already exists
        if (!index_to_region.emplace(cr.begin, cr).second) {
            THROW_ERROR("Duplicate starting index `" + std::to_string(cr.begin) + "` found in closed regions");
        }
        if (cr.begin < min_index) min_index = cr.begin;
        if (cr.begin > max_index) max_index = cr.begin;
    }

    std::vector<ClosedRegion> sorted_regions;
    sorted_regions.reserve(closed_regions.size());
    
    for (size_t i = min_index; i <= max_index; ++i){
        // checks if bucket exists, then stores the cr from bucket into sorted_regions
        std::unordered_map<size_t, ClosedRegion>::iterator it = index_to_region.find(i);
        if (it != index_to_region.end()) {
            sorted_regions.push_back(it->second);
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
              << "  unpaired=" << node->exclusive_unpaired_bases_count
              << "  children=" << node->children.size() << '\n';

    if (debug) {
        std::cout << *node << '\n';
    }

    for (const auto& child : node->children) {
        print_tree(child, depth + 1, debug);
    }
}

}  // namespace knotergy