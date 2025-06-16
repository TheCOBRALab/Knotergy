#include "Loop.hpp"

#include "../helpers/common.hpp"
#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
Loop::Loop(Bands& bands)
    : bands_(bands), entry(bands.entry_), regular_pairs(bands.entry_.get_regular_pairs()), pseudo_pairs(bands.entry_.get_pseudo_pairs())  {}

void Loop::build_tree(const std::vector<Pair>& pairs) {
    NULL;
}
}  // namespace compute_energy