#include "Loop.hpp"

#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"
#include "../helpers/common.hpp"

namespace compute_energy {
Loop::Loop(size_t start, size_t end, RNAEntry& rna_entry, Bands& bands,
           std::vector<Region>& stacks) {
    size_t start_ = start;
    size_t end_ = end;
    RNAEntry rna_entry_ = rna_entry;
    Bands bands_ = bands;
    std::vector<Region>& stacks_ = stacks;
}
}  // namespace ComputeEnergy