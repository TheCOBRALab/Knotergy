#pragma once

#include "../helpers/common.hpp"
#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {

class Loop {
   public:
    Loop(Bands& bands);
    ~Loop() = default;
    Bands& bands_;
    RNAEntry& entry;
    const std::vector<size_t>& pairings;

   private:
    void build_tree(const std::vector<Pair>& pairs);
    void add_loop(size_t i, size_t j);
};
}  // namespace compute_energy