#pragma once

#include <iostream>

#include "../loop_tree/LoopNode.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class PseudoknotFunctions {
   public:
    PseudoknotFunctions() = default;
    ~PseudoknotFunctions() = default;

    // Should be moved to private
    const int ext_pk_init_penalty = -138;                    // exterior pseudoloop initialization penalty
    const int pk_in_multi_penalty = 1007;                    // pseudoknot in multiloop penalty
    const int pk_in_pk_penalty = 1500;                       // pseudoknot in pseudoloop penalty
    const int band_penalty = 246;
    const int unpaired_in_pk_penalty = 6;                    // unpaired bases in pseudoknot penalty
    const int nested_cr_penalty = 96;                        // nested closed region penalty
    double pk_stack_penalty_x = 0.89;                        // stacked pair that spans a band penalty multiplier (band_energy * penalty)
    double pk_internal_penalty_x = 0.74;                     // internal pair that spans a band penalty multiplier (internal_energy * penalty)
    [[maybe_unused]] int pk_multi_init_penalty = 341;        // multiloop that spans a band penalty
    [[maybe_unused]] int pk_multi_bp_penalty = 56;           // base pair for multiloop that spans a band penalty
    [[maybe_unused]] int pk_unpaired_in_multi_penalty = 12;  // unpaired bases in a multiloop that spans a band penalty


    double pseudoknot_energy(const LoopNode& node, const std::string& sequence, RNAProcessedEntry processed_rna, bool round = false) {
        double energy = 0;
        int real_unpaired = node.exclusive_unpaired_bases_count;

        energy += init_penalty(node);
        for (Band band : node.bands){
            real_unpaired -= processed_rna.get_unpaired_count(band.left_border(), band.left_inner());
            real_unpaired -= processed_rna.get_unpaired_count(band.right_inner(), band.right_border());
        }

        energy += band_penalty * node.number_of_bands;
        std::cout << "Band penalty: "<< band_penalty * node.number_of_bands << std::endl;

        energy += unpaired_in_pk_penalty * real_unpaired;
        std::cout << "Unpaired penalty: "<< unpaired_in_pk_penalty * real_unpaired<< std::endl;

        energy += nested_cr_penalty * node.number_of_crossband_children;
        std::cout << "Nested penalty: "<< nested_cr_penalty * node.number_of_crossband_children << std::endl;

        energy += stack_and_internal_energy(node, sequence, round);

        for (std::shared_ptr<LoopNode> c : node.children){
            if (c->pseudo_type == PseudoNestedType::InsideBand){
                energy += pk_multi_bp_penalty * c->number_of_bands;
                std::cout << "multi bp penalty: "<< pk_multi_bp_penalty * c->number_of_bands << std::endl;
            }
        }

        // energy += pk_in_pk_penalty * (node.number_of_bands - 2);

        return energy;
    }

   private:
    ViennaFunctions vienna;

    double init_penalty(const LoopNode& node) {
        // initialization penalties
        double energy = 0;
        if (std::shared_ptr<LoopNode> parent = node.parent.lock()) {
            switch (parent->loop_type) {
                case (LoopType::External):
                    energy += ext_pk_init_penalty;
                    break;
                case (LoopType::Multibranch):
                    energy += pk_in_multi_penalty;
                    break;
                case (LoopType::Pseudoknot):
                    energy += node.pseudo_type == PseudoNestedType::InsideBand ? pk_in_multi_penalty : pk_in_pk_penalty; 
                    break;
                default:
                    std::cerr << "Warning: Parent of this node is not a pseudoknot, external, or multiloop" << node << std::endl;
                    break;
            }
        } else {
            THROW_ERROR("Parent node of pseudoknot (" + 
                        std::to_string(node.begin) + ", " +
                        std::to_string(node.end) + ") has expired.");
        }
        std::cout << "Init penalty: "<< energy << std::endl;
        return energy;
    }

    double stack_and_internal_energy(const LoopNode& node, const std::string& sequence, bool round) {
        double energy = 0;
        double offset = 0;

        for (const Band& band : node.bands) {
            const std::vector<BasePair>& bp = band.base_pairs();
            for (size_t idx = 0; idx < bp.size() - 1; ++idx) {

                // stack energy
                if (bp[idx].is_stack(bp[idx + 1])) {
                    int stack_energy = vienna.stack_energy(bp[idx], bp[idx + 1], sequence);
                    double stack_penalty = pk_stack_penalty_x * stack_energy;
                    energy += stack_penalty;
                    if (round) offset += stack_penalty - std::round(stack_penalty);
                    std::cout << "Stack penalty: "<< stack_penalty << std::endl;
                    continue;
                }

                // internal loop energy
                int internal_energy = vienna.internal_loop_energy(bp[idx], bp[idx + 1], sequence);
                double internal_penalty = pk_internal_penalty_x * internal_energy;
                energy += internal_penalty;
                if (round) offset += internal_penalty - std::round(internal_penalty);
                std::cout << "Internal penalty: "<< internal_penalty << std::endl;
            }
        }

        if (round){
            std::cout << "Offset: " << offset << std::endl;
            return energy - offset;
        }

        return energy;
        
    }
};

}  // namespace knotergy

// Test Cases

// CAGGGGAUAUUUUUCUUACUUAGCCAAACCUCCACCAACUCCGCCUGCUGGGCAACAAUCCUGAAGUGCGAGAGGCAUUAUAUUGAAUCCUGGUUCCAUAUUUCGACGAUAAAGCCAGGCUGGCGGACGGACCGACAGCAUUGAGAAACACACAUUGAAGUAGCGGUGGUUCGAAGACUUACGCUGAUUUGCGGGAGACGCACUGUUACUAUCACGUCCUGUUAUGGUUACUUAUUAGCCAGAUCAAGAC
// ..((((.......................))))......((((((.(((.(((....(((..(((.(((.....)))..((((.((((....)))).)))))))...)))...))).))).))))))[[..[[[.[[[[...[[[......[[.....((((((((((]].......]]]..]]]].....]]]....]])))))))))).....(((.((...((((((.....))))))...)).)))
// -39.19


// --------------------WORKING---------------------------

// GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC
// [[[[[.......((((((((((........]]]]]......))))))))))
// -20.16 HFold, -20.19 Knotergy
// HFold rounds every value, so it's off due to HFold's rounding error


// GGGGGGGAGGGGGAAAACCCCCAGGGGGGGGACCCCCCCAAACCCCCCCC
// (((((((.(((((....))))).[[[[[[[[.)))))))...]]]]]]]]
// hfold: -25.86, Knotergy: -25.89

// AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC
// .[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..
// ..........................((((((.......)))))).((((........................))))..

// HFold: -8.98, Knotergy: -8.98