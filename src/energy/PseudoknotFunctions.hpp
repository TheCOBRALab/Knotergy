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
        std::cout << "initial unpaired: " << unpaired << std::endl;
        for (Band band : node.bands) {
            unpaired -= processed_rna.get_unpaired_count(band.left_border(), band.left_inner());
            unpaired -= processed_rna.get_unpaired_count(band.right_inner(), band.right_border());

            // std::cout << "Band: (" << band.left_border() << ", " << band.right_border() << ")" <<std::endl;
            // std::cout << "Band: (" << band.left_inner() << ", " << band.right_inner() << ")" <<std::endl;
            // std::cout << "Remove: " << processed_rna.get_unpaired_count(band.left_border(), band.left_inner()) << std::endl;
            // std::cout << "Remove: " << processed_rna.get_unpaired_count(band.right_inner(), band.right_border()) << std::endl;
        }

        // these bases were removed twice (By using exclusive unpaired, and removing all base pairs in band)
        // so we're re-adding them
        for (std::shared_ptr<LoopNode> child : node.children){
            if (child->pseudo_type == PseudoNestedType::WithinBand){
                unpaired += child->total_unpaired_bases_count;
                // std::cout << "Add up: " << child->total_unpaired_bases_count << std::endl;
            }
        }

        double energy = 0;

        energy += init_penalty(node);

        energy += band_penalty * node.number_of_bands;
        std::cout << "Band penalty(" << node.begin << ", " << node.end << "): " << band_penalty * node.number_of_bands << std::endl;

        energy += unpaired_in_pk_penalty * unpaired;
        std::cout << "Unpaired penalty(" << unpaired << "): " << unpaired_in_pk_penalty * unpaired << std::endl;

        energy += nested_cr_penalty * node.number_of_nested_children;
        std::cout << "Nested penalty(" << node.begin << ", " << node.end << "): "<< nested_cr_penalty * node.number_of_nested_children
                  << std::endl;

        
        energy += loop_penalties(node, sequence, processed_rna, round);

        for (std::shared_ptr<LoopNode> c : node.children) {
            if (c->pseudo_type == PseudoNestedType::WithinBand) {
                energy += pk_multi_bp_penalty * c->number_of_bands;
                std::cout << "PKMloop bp penalty(" <<c->begin << ", " << c->end << "): " << pk_multi_bp_penalty * c->number_of_bands + 2 << std::endl;
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
                    energy += node.pseudo_type == PseudoNestedType::WithinBand ? pk_in_multi_penalty
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

        for (const Band& band : node.bands) {
            const std::vector<BasePair>& bps = band.base_pairs();
            const size_t n = bps.size();
            std::cout << band << std::endl;
            if (n < 2) THROW_ERROR("Less than 2 bands in pseudoknot");

            // loops through each base pair in band (except last one)
            for (size_t idx = 0; idx + 1 < n; ++idx) {
                const BasePair& bp = bps[idx];
                const BasePair& next_bp = bps[idx + 1];

                if (bp.is_stack(next_bp)) {
                    energy += pk_stack_energy(bp, next_bp, sequence, round);
                } else if (!bp.children.empty()) {
                    energy += pk_multiloop_energy(bp, next_bp, sequence, processed_rna);
                } else {
                    energy += pk_internal_energy(bp, next_bp, sequence, round);
                }
            }
        }
        return energy;
    }

    [[nodiscard]] double pk_stack_energy(const BasePair& bp, const BasePair& next_bp, const std::string& sequence, const bool& round){
        double stack_penalty = vienna.stack_energy(bp, next_bp, sequence) * pk_stack_penalty_x ;
        if (round) stack_penalty = std::round(stack_penalty);
        std::cout << "Stack penalty" << bp << ": "  << stack_penalty << std::endl;
        return stack_penalty;
    }

    [[nodiscard]] double pk_internal_energy(const BasePair& bp, const BasePair& next_bp, const std::string& sequence, const bool& round){
        double internal_penalty = vienna.internal_loop_energy(bp, next_bp, sequence) * pk_internal_penalty_x ;
        if (round) internal_penalty = std::round(internal_penalty);
        std::cout << "Internal penalty" << bp << ": " << internal_penalty << std::endl;
        return internal_penalty;
    }

     [[nodiscard]] double pk_multiloop_energy(const BasePair& bp, const BasePair& next_bp, [[maybe_unused]] const std::string& sequence, const ProcessedRNAEntry& processed_rna){
        double multiloop_penalty = pk_multi_init_penalty;
        std::cout << "PKMLoop Init penalty" << bp << ": " << pk_multi_init_penalty << std::endl;

        // for MLoop base pair, and nested basepair
        multiloop_penalty += pk_multi_bp_penalty * 2;
        std::cout << "PKMloop bp penalty" << bp << ": "<< pk_multi_bp_penalty * 2 << std::endl;

        // get unpaired count
        int unpaired = processed_rna.get_unpaired_count(bp.i, next_bp.i);
        unpaired += processed_rna.get_unpaired_count(next_bp.j, bp.j);
        for (BasePair child_bp : bp.children){
            unpaired -= processed_rna.get_unpaired_count(child_bp.i, child_bp.j);
        }

        // get unpaired penalty
        int pk_mloop_unpaired_energy = unpaired * pk_unpaired_in_multi_penalty;
        multiloop_penalty += pk_mloop_unpaired_energy;
        std::cout << "PKMLoop Unpaired penalty" << bp << ": " << pk_mloop_unpaired_energy << std::endl;

        // // get TerminalAU penalty
        // if ((sequence[bp.i] == 'A' && sequence[bp.j] == 'U') ||
        //     (sequence[bp.i] == 'U' && sequence[bp.j] == 'A')){
        //         multiloop_penalty += vienna.get_parameters()->TerminalAU;
        //         std::cout << "PKMloop Terminal AU Penalty: " << vienna.get_parameters()->TerminalAU << std::endl;
        // }

        return multiloop_penalty;
    }

};

}  // namespace knotergy

// Test Cases

// AAAAGGGAAAGGGUUUGGGAAAAGGGGGGGGGGGUUUUGGUUUUGGGGGGGGGGGGGGGGGAAAAGGGAAAACCCCCCCCCCCUUUUGGGGGAAAGGGUUUGGUUUU
// ((((xxx(((xxx)))xxx((((...........))))xx))))xxxxxxxxxxxxxxxxx((((xxx((((...........))))xxxxx(((xxx)))xx))))
// ((((...(((...)))...(((([[[[[[[[[[[))))..)))).................((((...((((]]]]]]]]]]])))).....(((...)))..)))) 
// -1.86


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


// 1 multiloop
// AAAAAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC
// (((((((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]
// 
// -0.74

// AAAGGAAAGGGUUUGGGGGGGGGGGAAAGGGUUUGGGGGGGGGGGUUUGGGGGCCCCCCCCCCCCCCCCCC
// (((..(((...)))...[[[[[...(((...)))...[[[[[...)))......]]]]]]]]]].......
// 0.7

// 2 multiloops
// AAAAGGAAAGGGGUUUGGGAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC
// ((((xx(((xxxx)))xxx(((xxx(((xxx)))xxx(((xx.......xxx)))xxx)))xxx))))xxx.......
// ((((..(((....)))...(((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]
// 4.65