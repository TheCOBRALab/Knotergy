#include "LoopFactory.hpp"

namespace compute_energy {

LoopFactory::LoopFactory(const RNAEntry& entry)
    : entry_(entry) {
        closed_regions_ = entry_.get_closed_regions();

        // Should already be sorted, this is just a safety pre-caution
        // Sorted where end index is in ascending order
        if (!std::is_sorted(closed_regions_.begin(), closed_regions_.end())) {
            std::sort(closed_regions_.begin(), closed_regions_.end());
        }

        build_tree();
    }

    void LoopFactory::build_tree() {
        std::stack<std::shared_ptr<LoopNode>> node_stack;

        root_node = std::make_shared<LoopNode>(ClosedRegion{NULL_INDEX, entry_.get_structure().size()});
        root_node -> loop_type = LoopType::External;
        node_stack.push(root_node);

        for (const ClosedRegion& closed_region: closed_regions_) {
            std::shared_ptr<LoopNode> node = std::make_shared<LoopNode>(closed_region);
            node -> num_of_unpaired_bases = entry_.get_unpaired_count(closed_region);
            
            if (node_stack.size() == 1){
                node_stack.push(node);
                node -> loop_type = get_loop_type(*node);
                // TODO: Find Bands
                // TODO: PseudoNested
                continue;
            }

            if ((node -> end > node_stack.top() -> end) || (node_stack.top()->begin == NULL_INDEX)){
                while (node_stack.size() > 1){
                    node->children.push_back(node_stack.top());
                    node_stack.pop();
                }
            }
            node_stack.push(node);
        }
        while (node_stack.size() > 1){
                    root_node->children.push_back(node_stack.top());
                    node_stack.pop();
                }
    }


void LoopFactory::print_tree() const {
    print_tree(root_node, 0);
}

void LoopFactory::print_tree(const std::shared_ptr<LoopNode>& node, int depth) const {
    std::cout << std::string(depth, '.') // indent with dots
        << '[' << node->begin
        << ',' << node->end
        << "]  "
        << "  unpaired=" << node->num_of_unpaired_bases
        << "  children=" << node->children.size()
        << '\n';

    for (const auto& child : node->children) {
        print_tree(child, depth + 1);
    }
}
    
    LoopType LoopFactory::get_loop_type(const LoopNode& node) {
        const std::vector<size_t>& pairings = entry_.get_pairings();
        
        if (pairings[node.begin] != node.end){
            return LoopType::Pseudoknot;
        }

        // Zero Children
        if (node.children.empty()){
            return LoopType::Hairpin;
        }

        // if only one child and it's not pseudoknotted...
        if (node.children.size() == 1 && node.children[0]->loop_type != LoopType::Pseudoknot){
            // ... and if child's base pair is stacked with this loop's pair, then stacked
            if ((node.children[0]->begin == node.begin + 1) && node.children[0]->end == node.end -1){
                return LoopType::Stack;
            }
            // ...and if not stacked, then bulge or internal loop
            return LoopType::Internal;
        }

        // Multiple children or a pseudoknotted child = Multiloop
        return LoopType::Multi;

    }


}