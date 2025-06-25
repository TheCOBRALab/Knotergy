#pragma once

#include <cstddef>
#include <vector>

#include "../helpers/common.hpp"
#include "../loops/LoopNode.hpp"
#include "Band.hpp"
namespace compute_energy {

// walk one step along a perfect stack:  (i+1) pairs (j−1)
class BandUtils {
   public:
    BandUtils(const std::vector<size_t>& pairings)
        : pairings_(pairings), done_(pairings_.size(), false) {};

    std::vector<Band> find_bands_in_region(size_t left, size_t right);
    std::vector<Band> find_bands_in_region() {
        return find_bands_in_region(0, pairings_.size() - 1);
    }
    void annotate_bands(const std::shared_ptr<LoopNode>& node);

   private:
    const std::vector<size_t> pairings_;
    std::vector<bool> done_;
    bool extend_stem(size_t& il, size_t& jr);
};

}  // namespace compute_energy