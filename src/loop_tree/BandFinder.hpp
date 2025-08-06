#pragma once
#include <iostream>

#include "LoopNode.hpp"
#include "../preprocessing/RNAProcessedEntry.hpp"
namespace knotergy {
class BandFinder {
   public:
    BandFinder();
    static std::vector<Band> find_bands(const size_t& left_bound, 
                                        const size_t& right_bound,
                                        const LoopType& loop_type,
                                        const std::vector<size_t>& pairings,
                                        const std::vector<size_t>& cr_pairings);

    static std::vector<Band> find_bands(const LoopNode& node, const RNAProcessedEntry& processed_rna);

   private:
    static bool extend_stem(size_t& i_prime,
                            size_t& j_prime,
                            const size_t& left_bound,
                            const size_t& right_bound,
                            const std::vector<size_t>& pairings);

    static std::unordered_map<size_t, size_t> const generate_next_bp_map(const size_t& left_bound,
                                                                   const size_t& right_bound,
                                                                   const std::vector<size_t>& pairings,
                                                                   const std::vector<size_t>& cr_pairings);

    static std::unordered_map<size_t, size_t> const generate_next_bp_map(const LoopNode& node, const RNAProcessedEntry& processed_entry);

    static std::unordered_map<size_t, size_t> const generate_prev_bp_map(const size_t& left_bound,
                                                                   const size_t& right_bound,
                                                                   const std::vector<size_t>& pairings,
                                                                   const std::vector<size_t>& cr_pairings);

    static std::unordered_map<size_t, size_t> const generate_prev_bp_map(const LoopNode& node, const RNAProcessedEntry& processed_entry);
};
}  // namespace knotergy