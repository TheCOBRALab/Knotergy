#pragma once

#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"
#include "../helpers/common.hpp"

namespace compute_energy {
class Loop {
   public:
    Loop(size_t start, size_t end, RNAEntry& rna_entry, Bands& bands, std::vector<Region>& stacks);
    ~Loop();
};
}  // namespace ComputeEnergy