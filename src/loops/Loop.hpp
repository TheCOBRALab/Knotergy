#pragma once

#include "Bands.hpp"
#include "ComputeEnergy.hpp"
#include "RNAEntry.hpp"

namespace ComputeEnergy {
class Loop {
   public:
    Loop(size_t start, size_t end, RNAEntry& rna_entry, Bands& bands, std::vector<Region>& stacks);
    ~Loop();
};
}  // namespace ComputeEnergy