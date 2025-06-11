#include "Loop.hpp"

#include "../helpers/common.hpp"
#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
Loop::Loop(Bands& bands,
           std::vector<Region>& stacks) {
    bands_ = bands;
    stacks_ = stacks;
    entry_ = bands_.entry_;
}

void Loop::build_tree(){
    const std::vector<size_t> pairings = entry_.get_pairings();
    size_t left, right;
    for (size_t i = 0; i <= pairings.size(); ++i){
        if (pairings[i] != NULL_INDEX){
            NULL; // TODO
        }
    }
}

bool Loop::Add(size_t a, int&b, int &e){
    const std::vector<size_t> pairings = entry_.get_pairings();

    if (pairings[a] == NULL_INDEX){  // a is unpaired; do nothing
		return false;
    };

    if (a < pairings[a]) { // potentially closed region [a,bp(a)] added to stack
        Region current{a, pairings[a]};
        stacks_.push_back(current);
        return false;
    }

    if (pairings.back() == pairings[a])

}
}  // namespace compute_energy