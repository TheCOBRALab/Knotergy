#pragma once
#include <iostream>
#include <utility>

#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "LoopNode.hpp"
namespace knotergy {

struct BandLink {
    size_t value = NULL_INDEX;
    size_t prev = NULL_INDEX;
    size_t next = NULL_INDEX;
};

class BandFinder {
   public:
    BandFinder();
    static std::vector<Band> find_bands(const size_t& left_bound, const size_t& right_bound,
                                        const LoopType& loop_type,
                                        const std::vector<size_t>& pairings,
                                        const std::vector<size_t>& cr_pairings);

    static std::vector<Band> find_bands(const LoopNode& node,
                                        const ProcessedRNAEntry& processed_rna);

   private:
    static bool extend_stem(size_t& i_prime, size_t& j_prime,
                            const std::unordered_map<size_t, BandLink>& aux_bands,
                            const std::vector<size_t>& pairings);

    static std::unordered_map<size_t, BandLink> const generate_band_links(
        const size_t& left_bound, const size_t& right_bound, const std::vector<size_t>& pairings,
        const std::vector<size_t>& cr_pairings);

    static std::unordered_map<size_t, BandLink> const generate_band_links(
        const LoopNode& node, const ProcessedRNAEntry& processed_entry);
};
}  // namespace knotergy