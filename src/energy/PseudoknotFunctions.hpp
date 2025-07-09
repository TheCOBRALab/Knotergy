#pragma once

#include <iostream>

#include "../loop_tree/LoopNode.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class PseudoknotFunctions {
   public:
    PseudoknotFunctions() = default;
    ~PseudoknotFunctions() = default;
    // GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC
    // [[[[[.......((((((((((........]]]]]......))))))))))
    // -20.16
    double pseudoknot_energy(const LoopNode& node, const std::string& sequence) {

        const int ext_pk_init_penalty = -138;  // exterior pseudoloop initialization penalty
        [[maybe_unused]] const int pk_in_multi_penalty = 1007; // pseudoknot in multiloop penalty
        const int pk_in_pk_penalty = 1500;    // pseudoknot in pseudoloop penalty
        const int band_penalty = 246;
        const int unpaired_in_pk_penalty = 6; // unpaired bases in pseudoknot penalty
        const int nested_cr_penalty = 96; // nested closed region penalty 
        
        double pk_stack_penalty_x = 0.89; // stacked pair that spans a band penalty multiplier (band_energy * penalty)
        double pk_internal_penalty_x = 0.74; // internal pair that spans a band penalty multiplier (internal_energy * penalty)

        [[maybe_unused]] int multi_init_penalty = 339; // multiloop initialization penalty
        [[maybe_unused]] int multi_bp_penalty = 3; // multiloop base pair penalty
        [[maybe_unused]] int unpaired_in_multi_penalty = 2; // unpaired base in multiloop penalty

        [[maybe_unused]] int pk_multi_init_penalty = 341; // multiloop that spans a band penalty
        [[maybe_unused]] int pk_multi_bp_penalty = 56; //base pair for multiloop that spans a band penalty
        [[maybe_unused]] int pk_unpaired_in_multi_penalty = 12; // unpaired bases in a multiloop that spans a band penalty


        double energy = 0;

        // std::cout << node << std::endl;
        energy += ext_pk_init_penalty;

        energy += pk_in_pk_penalty * (node.number_of_bands - 2);
        std::cout << "pk_in_pk: " << pk_in_pk_penalty * (node.number_of_bands - 2) << std::endl;

        energy += band_penalty * node.number_of_bands;
        std::cout << "band_penalty: " << band_penalty * node.number_of_bands << std::endl;

        energy += unpaired_in_pk_penalty * node.number_of_exclusive_unpaired_bases;
        std::cout << "unpaired_in_pk_penalty: " << unpaired_in_pk_penalty * node.number_of_exclusive_unpaired_bases << std::endl;

        energy += nested_cr_penalty * node.number_of_children_outside_band;
        std::cout << "nested_cr_penalty: " << nested_cr_penalty * node.number_of_children_outside_band << std::endl;


        for (const Band& band : node.bands) {
            const std::vector<Pair>& bp = band.base_pairs();
            for (size_t idx = 0; idx < bp.size() - 1; ++idx) {
                if (bp[idx].is_stack(bp[idx+1])){
                    energy += pk_stack_penalty_x * vienna.stack_energy(bp[idx], bp[idx + 1], sequence);
                    continue;
                }
                energy += pk_internal_penalty_x * vienna.internal_loop_energy(bp[idx], bp[idx + 1], sequence);
            }
        }

        return energy;
    }

   private:
    ViennaFunctions vienna;
};
}  // namespace knotergy