#pragma once

#include <vector>

#include "ClosedRegion.hpp"
#include "RNAEntry.hpp"
#include "ProcessedRNAEntry.hpp"
namespace knotergy {
class RNAProcessor {
   public:
    RNAProcessor();
    static ProcessedRNAEntry process_rna(RNAEntry rna);

   private:
    static std::vector<size_t> compute_pairings(RNAEntry rna);
    static std::vector<ClosedRegion> compute_closed_regions(std::vector<size_t> pairings);
    static std::vector<size_t> compute_closed_regions_pairings(std::vector<ClosedRegion> closed_regions, size_t rna_size);
    static std::vector<int> compute_unpaired_counts(std::vector<size_t> pairings);
};
}  // namespace knotergy