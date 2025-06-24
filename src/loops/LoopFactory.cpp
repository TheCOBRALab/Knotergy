#include "LoopFactory.hpp"

#include <algorithm>
#include <memory>
#include <stack>

namespace compute_energy {

LoopFactory::LoopFactory(const RNAEntry& entry) : entry_(entry) {
    closed_regions_ = entry_.get_closed_regions();

    // Should already be sorted, this is just a safety pre-caution
    // Sorted where end index is in ascending order
    [[likely]] if (!std::is_sorted(closed_regions_.begin(), closed_regions_.end())) {
        std::cerr << "Warning: Closed regions are not sorted by end index. Sorting now.\n";
        std::sort(closed_regions_.begin(), closed_regions_.end());
    }

    build_tree();
}

void LoopFactory::build_tree() {
    std::stack<std::shared_ptr<LoopNode>> node_stack;

    root_node_ =
        std::make_shared<LoopNode>(ClosedRegion{NULL_INDEX, entry_.get_structure().size()});
    root_node_->loop_type = LoopType::External;
    node_stack.push(root_node_);

    for (const ClosedRegion& closed_region : closed_regions_) {
        std::shared_ptr<LoopNode> node = std::make_shared<LoopNode>(closed_region);
        node->num_of_unpaired_bases = entry_.get_unpaired_count(closed_region);

        if (node_stack.size() == 1) {
            node_stack.push(node);
            node->loop_type = get_loop_type(*node);
            // TODO: Find Bands
            // TODO: PseudoNested
            continue;
        }

        // Pop nodes from the stack until we find a node that is not a child of the current node
        if ((node->begin < node_stack.top()->begin) || (node_stack.top()->begin == NULL_INDEX)) {
            while (node_stack.size() > 1) {
                node->children.push_back(node_stack.top());
                node_stack.pop();
            }
        }
        node_stack.push(node);
    }

    // Pop all remaining nodes in the stack and add them to the root node
    while (node_stack.size() > 1) {
        root_node_->children.push_back(node_stack.top());
        node_stack.pop();
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

void LoopFactory::print_tree() const {
    std::cout << "Root Node: ";
    print_tree(root_node_, 0);
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