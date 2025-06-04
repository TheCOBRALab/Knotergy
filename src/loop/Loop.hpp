#pragma once

#include "../helpers/common.hpp"
#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
class Loop {
   public:
    Loop(size_t start, size_t end, RNAEntry& rna_entry, Bands& bands, std::vector<Region>& stacks);
    ~Loop();
};
}  // namespace compute_energy