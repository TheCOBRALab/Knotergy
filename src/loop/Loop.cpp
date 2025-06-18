#include "Loop.hpp"

#include "../helpers/common.hpp"
#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
Loop::Loop(Bands& bands) : bands_(bands), entry(bands.entry_) {}

void Loop::build_tree(const std::vector<Region>& pairs) {
    NULL;
}
}  // namespace compute_energy