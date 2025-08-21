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

    // multiloop that spans a band penalty
    int pk_multi_init_penalty = 341;
    // base pair for multiloop that spans a band penalty
    int pk_multi_bp_penalty = 56;
    // unpaired bases in a multiloop that spans a band penalty
    int pk_unpaired_in_multi_penalty = 12;

    double pseudoknot_energy(const LoopNode& node, const std::string& sequence,
                             ProcessedRNAEntry processed_rna, bool round = false) {
        
        // Unpaired within bands are already included in stack_and_internal_energy
        int unpaired = node.exclusive_unpaired_bases_count;
        for (Band band : node.bands) {
            unpaired -= processed_rna.get_unpaired_count(band.left_border(), band.left_inner());
            unpaired -= processed_rna.get_unpaired_count(band.right_inner(), band.right_border());
        }

        double energy = 0;

        energy += init_penalty(node);

        energy += band_penalty * node.number_of_bands;
        std::cout << "Band penalty: " << band_penalty * node.number_of_bands << std::endl;

        energy += unpaired_in_pk_penalty * unpaired;
        std::cout << "Unpaired penalty: " << unpaired_in_pk_penalty * unpaired << std::endl;

        energy += nested_cr_penalty * node.number_of_crossband_children;
        std::cout << "Nested penalty: " << nested_cr_penalty * node.number_of_crossband_children
                  << std::endl;

        
        energy += loop_penalties(node, sequence, processed_rna, round);

        for (std::shared_ptr<LoopNode> c : node.children) {
            if (c->pseudo_type == PseudoNestedType::InsideBand) {
                energy += pk_multi_bp_penalty * c->number_of_bands;
                std::cout << "PKMloop bp penalty: " << pk_multi_bp_penalty * c->number_of_bands
                          << std::endl;
            }
        }

        // energy += pk_in_pk_penalty * (node.number_of_bands - 2);

        return energy;
    }

   private:
    ViennaFunctions vienna;

    [[nodiscard]] double init_penalty(const LoopNode& node) {
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
                    energy += node.pseudo_type == PseudoNestedType::InsideBand ? pk_in_multi_penalty
                                                                               : pk_in_pk_penalty;
                    break;
                default:
                    std::cerr << "Warning: Parent of this node is not a pseudoknot, external, or "
                                 "multiloop"
                              << node << std::endl;
                    break;
            }
        } else {
            THROW_ERROR("Parent node of pseudoknot (" + std::to_string(node.begin) + ", " +
                        std::to_string(node.end) + ") has expired.");
        }
        std::cout << "Init penalty: " << energy << std::endl;
        return energy;
    }

    [[nodiscard]] double loop_penalties(const LoopNode& node, const std::string& sequence, const ProcessedRNAEntry& processed_rna, bool round) {
        double energy = 0;
        double offset = 0;

        for (const Band& band : node.bands) {
            const std::vector<BasePair>& bps = band.base_pairs();
            const size_t n = bps.size();
            if (n < 2) THROW_ERROR("Less than 2 bands in pseudoknot");

            for (size_t idx = 0; idx + 1 < n; ++idx) {
                const BasePair& bp = bps[idx];
                const BasePair& next_bp = bps[idx + 1];

                // stack energy
                if (bp.is_stack(next_bp)) {
                    energy += pk_stack_energy(bp, next_bp, sequence, round, offset);
                    continue;
                }

                if (!bp.children.empty()){
                    energy += pk_multi_init_penalty;
                    std::cout << "PKMLoop Init penalty: " << pk_multi_init_penalty << std::endl;

                    // get unpaired count
                    int unpaired = processed_rna.get_unpaired_count(bp.i, next_bp.i);
                    unpaired += processed_rna.get_unpaired_count(next_bp.j, bp.j);
                    for (BasePair child_bp : bp.children){
                        unpaired -= processed_rna.get_unpaired_count(child_bp.i, child_bp.j);
                    }

                    // get unpaired penalty
                    int pk_mloop_unpaired_energy = unpaired * pk_unpaired_in_multi_penalty;
                    energy += pk_mloop_unpaired_energy;
                    std::cout << "PKMLoop Unpaired penalty: " << pk_mloop_unpaired_energy << std::endl;

                    // get TerminalAU penalty
                    if ((sequence[bp.i] == 'A' && sequence[bp.j] == 'U') ||
                        (sequence[bp.i] == 'U' && sequence[bp.j] == 'A')){
                            energy += vienna.get_parameters()->TerminalAU;
                            std::cout << "PKMloop Terminal AU Penalty: " << vienna.get_parameters()->TerminalAU << std::endl;
                    }
                    continue;
                }

                // internal loop energy
                energy += pk_internal_energy(bp, next_bp, sequence, round, offset);
            }
        }
     
        if (round) std::cout << "Offset: " << offset << std::endl;
        return round ? (energy - offset) : energy;
    }

    [[nodiscard]] double pk_stack_energy(const BasePair& bp, const BasePair& next_bp, const std::string& sequence, const bool& round, double& offset){
        int stack_energy = vienna.stack_energy(bp, next_bp, sequence);
        double stack_penalty = pk_stack_penalty_x * stack_energy;
        if (round) offset += stack_penalty - std::round(stack_penalty);
        std::cout << "Stack penalty: " << stack_penalty << std::endl;
        return stack_penalty;
    }

    [[nodiscard]] double pk_internal_energy(const BasePair& bp, const BasePair& next_bp, const std::string& sequence, const bool& round, double& offset){
        int internal_energy = vienna.internal_loop_energy(bp, next_bp, sequence);
        double internal_penalty = pk_internal_penalty_x * internal_energy;
        if (round) offset += internal_penalty - std::round(internal_penalty);
        std::cout << "Internal penalty: " << internal_penalty << std::endl;
        return internal_penalty;
    }

};

}  // namespace knotergy

// Test Cases

// CAGGGGAUAUUUUUCUUACUUAGCCAAACCUCCACCAACUCCGCCUGCUGGGCAACAAUCCUGAAGUGCGAGAGGCAUUAUAUUGAAUCCUGGUUCCAUAUUUCGACGAUAAAGCCAGGCUGGCGGACGGACCGACAGCAUUGAGAAACACACAUUGAAGUAGCGGUGGUUCGAAGACUUACGCUGAUUUGCGGGAGACGCACUGUUACUAUCACGUCCUGUUAUGGUUACUUAUUAGCCAGAUCAAGAC
// ..((((.......................))))......((((((.(((.(((....(((..(((.(((.....)))..((((.((((....)))).)))))))...)))...))).))).))))))[[..[[[.[[[[...[[[......[[.....((((((((((]].......]]]..]]]].....]]]....]])))))))))).....(((.((...((((((.....))))))...)).)))
// -39.19

// 1 multiloop
// AAAAAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC
// (((((((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]
// 
// -0.5

// AAAGGAAAGGGUUUGGGGGGGGGGGAAAGGGUUUGGGGGGGGGGGUUUGGGGGCCCCCCCCCCCCCCCCCC
// (((..(((...)))...[[[[[...(((...)))...[[[[[...)))......]]]]]]]]]].......
// 0.13

// 2 multiloops
// AAAAGGAAAGGGGUUUGGGAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC
// ((((..(((....)))...(((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]
// 5.13

// AAAAGGGAAAGGGUUUGGGAAAAGGGGGGGGGGGUUUUGGUUUUGGGGGGGGGGGGGGGGGAAAAGGGAAAACCCCCCCCCCCUUUUGGGGGAAAGGGUUUGGUUUU
// ((((...(((...)))...(((([[[[[[[[[[[))))..)))).................((((...((((]]]]]]]]]]])))).....(((...)))..)))) 
// -1.38


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