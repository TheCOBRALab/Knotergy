#pragma once
#include "LoopNode.hpp"
#include <iostream>
namespace knotergy{
    class BandFinder{
        public:
        BandFinder();
        static std::vector<Band> find_bands(const size_t& left_bound, const size_t& right_bound, const LoopType& loop_type, const std::vector<size_t>& pairings, const std::vector<size_t>& cr_pairings);
        static std::vector<Band> find_bands(const LoopNode& node, const std::vector<size_t>& pairings, const std::vector<size_t>& cr_pairings);

        private:
        static bool extend_stem(size_t& i_prime, size_t& j_prime, const size_t& left_bound, const size_t& right_bound, const std::vector<size_t>& pairings);
    };
}