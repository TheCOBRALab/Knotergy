#include "LoopFactory.hpp"

#include <algorithm>
#include <memory>
#include <stack>

namespace compute_energy {

LoopFactory::LoopFactory(const RNAEntry& entry) : entry_(entry) {
    closed_regions_ = entry_.get_closed_regions();
    structure_length_ = entry_.get_structure().size();

    // Sorted by begin index in ascending order (Complexity: O(n))
    closed_regions_ = closed_region_bucket_sort(closed_regions_, structure_length_);

    build_tree();
}

void LoopFactory::build_tree() {
    std::stack<std::shared_ptr<LoopNode>> node_stack;

    root_node_ = std::make_shared<LoopNode>(ClosedRegion{NULL_INDEX, structure_length_});
    root_node_->loop_type = LoopType::External;
    node_stack.push(root_node_);

    for (const ClosedRegion& closed_region : closed_regions_) {
        // pop until node_stack.end() is parent of current node
        while (node_stack.top()->end < closed_region.begin) {
            node_stack.pop();
        }

        std::shared_ptr<LoopNode>& parent = node_stack.top();
        parent->children.emplace_back(std::make_shared<LoopNode>(closed_region));
        std::shared_ptr<LoopNode>& child = parent->children.back();

        child->num_of_unpaired_bases = entry_.get_unpaired_count(closed_region);
        child->loop_type = get_loop_type(*child);

        node_stack.push(child);
    }
}

LoopType LoopFactory::get_loop_type(const LoopNode& node) {
    const std::vector<size_t>& pairings = entry_.get_pairings();

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
    return LoopType::Multi;
}

void LoopFactory::PseudoNestedCheck(const LoopNode& node) {
    if (node.loop_type != LoopType::Pseudoknot) return;

    for (std::shared_ptr<LoopNode> n : node.children) {
    }
}

std::vector<ClosedRegion> LoopFactory::closed_region_bucket_sort(
    std::vector<ClosedRegion>& closed_regions, size_t structure_length) {
    /*  one linear pass buckets the regions by their left index  */
    std::vector<std::vector<ClosedRegion>> buckets(structure_length);
    for (const ClosedRegion& cr : closed_regions) {
        if (!buckets[cr.begin].empty()) {
            throw std::runtime_error("Duplicate starting index " + std::to_string(cr.begin) +
                                     "found in closed regions");
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

void LoopFactory::print_tree() const {
    for (const auto& child : root_node_->children) {
        print_tree(child, 0);
    }
}

void LoopFactory::print_tree(const std::shared_ptr<LoopNode>& node, size_t depth) const {
    std::cout << std::string(depth, '.')  // indent with dots
              << '[' << node->begin << ',' << node->end << "]  "
              << "  unpaired=" << node->num_of_unpaired_bases
              << "  children=" << node->children.size() << '\n';

    for (const auto& child : node->children) {
        print_tree(child, depth + 1);
    }
}

}  // namespace compute_energy