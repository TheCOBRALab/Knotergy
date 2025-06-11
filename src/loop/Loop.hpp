#pragma once

#include "../helpers/common.hpp"
#include "../rna_regions/Bands.hpp"
#include "../rna_regions/RNAEntry.hpp"

namespace compute_energy {
class Loop {
   public:
    Loop(Bands& bands, std::vector<Region>& stacks);
    ~Loop();
    Bands bands_;
    std::vector<Region> stacks_;
    RNAEntry entry_;

    private:
    void build_tree();
    bool Add(size_t a, int&b, int &e);
};
}  // namespace compute_energy