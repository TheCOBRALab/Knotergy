#pragma once

#include <iostream>

#include "../loop_tree/LoopNode.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class PseudoknotFunctions {
   public:
    PseudoknotFunctions() = default;
    ~PseudoknotFunctions() = default;

    double pseudoknot_energy(const LoopNode& node, const std::string& sequence) {
        const int ext_pk_init_penalty = -138;  // exterior pseudoloop initialization penalty
        const int pk_in_multi_penalty = 1007;  // pseudoknot in multiloop penalty
        const int pk_in_pk_penalty = 1500;     // pseudoknot in pseudoloop penalty
        const int band_penalty = 246;
        const int unpaired_in_pk_penalty = 6;  // unpaired bases in pseudoknot penalty
        const int nested_cr_penalty = 96;      // nested closed region penalty

        // stacked pair that spans a band penalty multiplier (band_energy * penalty)
        double pk_stack_penalty_x = 0.89;
        // internal pair that spans a band penalty multiplier (internal_energy * penalty)
        double pk_internal_penalty_x = 0.74;

        [[maybe_unused]] int pk_multi_init_penalty = 341;  // multiloop that spans a band penalty
        [[maybe_unused]] int pk_multi_bp_penalty =
            56;  // base pair for multiloop that spans a band penalty
        [[maybe_unused]] int pk_unpaired_in_multi_penalty =
            12;  // unpaired bases in a multiloop that spans a band penalty

        double energy = 0;

        // std::cout << node << std::endl;
        energy += ext_pk_init_penalty;

        if (std::shared_ptr<LoopNode> parent = node.parent.lock()) {
            switch (parent->loop_type){
                case(LoopType::Multibranch):
                    energy += pk_in_multi_penalty;
                    break;
                case(LoopType::Pseudoknot):
                    energy += node.pseudo_type == PseudoNestedType::InsideBand ? pk_in_pk_penalty : pk_multi_bp_penalty;
                    energy += pk_in_multi_penalty;
                    break;
                default:
                    break;
            }
        } else {
            THROW_ERROR("Parent node of pseudoknot (" + std::to_string(node.begin) + ", " +
                        std::to_string(node.end) + ") has expired.");
        }

        // energy += pk_in_pk_penalty * (node.number_of_bands - 2);
        // std::cout << "pk_in_pk: " << pk_in_pk_penalty * (node.number_of_bands - 2) << std::endl;

        energy += band_penalty * node.number_of_bands;
        energy += unpaired_in_pk_penalty * node.number_of_exclusive_unpaired_bases;
        energy += nested_cr_penalty * node.number_of_crossband_children;

        for (const Band& band : node.bands) {
            const std::vector<BasePair>& bp = band.base_pairs();
            for (size_t idx = 0; idx < bp.size() - 1; ++idx) {
                if (bp[idx].is_stack(bp[idx + 1])) {
                    energy += pk_stack_penalty_x * vienna.stack_energy(bp[idx], bp[idx + 1], sequence);
                    continue;
                }
                energy += pk_internal_penalty_x *
                          vienna.internal_loop_energy(bp[idx], bp[idx + 1], sequence);
            }
        }

        return energy;
    }

   private:
    ViennaFunctions vienna;
};
}  // namespace knotergy

// Test Cases

// UUAAAAGGGAUGCCUCUCCUGUUCAUCUUGUGGAGAAGCAUUCGAUAAGGUCAUCAUAAUGGGUCCAGCUUUGCGACCUGGCGAGAUUAGUCAGGAAAAUGUGAAGUGGGUCUUCGCUUUCCA
// .......((((((.((((([[...[[[[[[[))))).)))))).]]]]]]]...]]....(((.((.(((((((..((((((.......)))))).....))))))).)).))).........
// -27.99

// AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC
// .[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..
// -8.98

// GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC
// [[[[[.......((((((((((........]]]]]......))))))))))
// -20.16