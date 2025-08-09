#pragma once
#include <iostream>
#include <utility>

#include "LoopNode.hpp"
#include "../preprocessing/RNAProcessedEntry.hpp"
namespace knotergy {

struct Aux_Band {
    size_t value = NULL_INDEX;
    size_t prev = NULL_INDEX;
    size_t next = NULL_INDEX;
    
    Aux_Band() = default;  
    Aux_Band (size_t v) : value{v} {}
    Aux_Band (size_t v, size_t p) : value{v}, prev{p} {}
    Aux_Band (size_t v, size_t p, size_t n) : value{v}, prev{p}, next{n} {}
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
                            const std::unordered_map<size_t, Aux_Band>& aux_bands,
                            const std::vector<size_t>& pairings);

    static std::unordered_map<size_t, Aux_Band> const generate_Aux_Band_list(const size_t& left_bound,
                                                                   const size_t& right_bound,
                                                                   const std::vector<size_t>& pairings,
                                                                   const std::vector<size_t>& cr_pairings);

    static std::unordered_map<size_t, Aux_Band> const generate_Aux_Band_list(const LoopNode& node, const RNAProcessedEntry& processed_entry);
    };
}  // namespace knotergy