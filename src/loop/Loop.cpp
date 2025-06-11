#include "Loop.hpp"

#include "../helpers/common.hpp"
#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
Loop::Loop(Bands& bands)
    : bands_(bands), entry(bands.entry_), pairings(bands.entry_.get_pairings()) {}

void Loop::build_tree(const std::vector<Pair>& pairs) {
    NULL;
}
}  // namespace compute_energy