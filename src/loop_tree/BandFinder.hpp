#pragma once
#include <iostream>
#include <utility>

#include "LoopNode.hpp"
#include "../preprocessing/RNAProcessedEntry.hpp"
namespace knotergy {

struct LinkedList {
    size_t value = NULL_INDEX;
    std::shared_ptr<LinkedList> next;
    std::weak_ptr<LinkedList> prev;

    LinkedList(size_t v): value{v} {};
};

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

    static std::pair<std::shared_ptr<LinkedList>, std::shared_ptr<LinkedList>> const generate_aux_band_list(const size_t& left_bound,
                                                                   const size_t& right_bound,
                                                                   const std::vector<size_t>& pairings,
                                                                   const std::vector<size_t>& cr_pairings);

    static std::pair<std::shared_ptr<LinkedList>, std::shared_ptr<LinkedList>> const generate_aux_band_list(const LoopNode& node, const RNAProcessedEntry& processed_entry);
    };
}  // namespace knotergy