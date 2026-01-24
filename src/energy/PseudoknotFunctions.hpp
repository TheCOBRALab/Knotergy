#pragma once

#include <iostream>

#include "../io/PseudoknotParams.hpp"
#include "../io/ViennaParams.hpp"
#include "../loop_tree/LoopNode.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class PseudoknotFunctions {
   public:
    static double pseudoknot_energy(const LoopNode& node, const std::string& sequence,
                                    ProcessedRNAEntry processed_rna, bool round = false) {
        // Unpaired within bands are already included in stack_and_internal_energy
        int unpaired = node.exclusive_unpaired_bases_count;
        // std::cout << "initial unpaired: " << unpaired << std::endl;
        for (Band band : node.bands) {
            unpaired -= processed_rna.get_unpaired_count(band.left_border(), band.left_inner());
            unpaired -= processed_rna.get_unpaired_count(band.right_inner(), band.right_border());

            // std::cout << "Band: (" << band.left_border() << ", " << band.right_border() << ")"
            // <<std::endl; std::cout << "Band: (" << band.left_inner() << ", " <<
            // band.right_inner() << ")" <<std::endl; std::cout << "Remove: " <<
            // processed_rna.get_unpaired_count(band.left_border(), band.left_inner()) << std::endl;
            // std::cout << "Remove: " << processed_rna.get_unpaired_count(band.right_inner(),
            // band.right_border()) << std::endl;
        }

        // these bases were removed twice (By using exclusive unpaired, and removing all base pairs
        // in band) so we're re-adding them
        for (std::shared_ptr<LoopNode> child : node.children) {
            if (child->pseudo_type == PseudoNestedType::WithinBand) {
                unpaired += child->total_unpaired_bases_count;
            }
        }

        double energy = 0;

        energy += PseudoknotFunctions::init_penalty(node);
        energy += PseudoknotParams::pkp->band * node.number_of_bands;
        energy += PseudoknotParams::pkp->unpaired_in_pk * unpaired;
        energy += PseudoknotParams::pkp->cr_in_pk * node.number_of_nested_children;
        energy += PseudoknotFunctions::loop_penalties(node, sequence, processed_rna, round);

        for (std::shared_ptr<LoopNode> c : node.children) {
            if (c->pseudo_type == PseudoNestedType::WithinBand) {
                energy += PseudoknotParams::pkp->pk_mloop_bp * c->number_of_bands;
            }
        }
        return energy;
    }

   private:
    [[nodiscard]] static double init_penalty(const LoopNode& node) {
        // initialization penalties
        double energy = 0;
        if (std::shared_ptr<LoopNode> parent = node.parent.lock()) {
            switch (parent->loop_type) {
                case (LoopType::External):
                    energy += PseudoknotParams::pkp->pk_in_ext;
                    break;
                case (LoopType::Multibranch):
                    energy += PseudoknotParams::pkp->pk_in_mloop;
                    break;
                case (LoopType::Pseudoknot):
                    energy += node.pseudo_type == PseudoNestedType::WithinBand
                                  ? PseudoknotParams::pkp->pk_in_mloop
                                  : PseudoknotParams::pkp->pk_in_pk;
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
        return energy;
    }

    [[nodiscard]] static double loop_penalties(const LoopNode& node, const std::string& sequence,
                                               const ProcessedRNAEntry& processed_rna, bool round) {
        double energy = 0;

        for (const Band& band : node.bands) {
            const std::vector<BasePair>& bps = band.base_pairs();
            const size_t n = bps.size();

            // loops through each base pair in band (except last one)
            for (size_t idx = 0; idx + 1 < n; ++idx) {
                const BasePair& bp = bps[idx];
                const BasePair& next_bp = bps[idx + 1];

                if (bp.is_stack(next_bp)) {
                    energy += PseudoknotFunctions::pk_stack_energy(bp, next_bp, sequence, round);
                } else if (!bp.children.empty()) {
                    energy += PseudoknotFunctions::pk_multiloop_energy(bp, next_bp, processed_rna);
                } else {
                    energy += PseudoknotFunctions::pk_internal_energy(bp, next_bp, sequence, round);
                }
            }
        }
        return energy;
    }

    [[nodiscard]] static double pk_stack_energy(const BasePair& bp, const BasePair& next_bp,
                                                const std::string& sequence, const bool& round) {
        double stack_penalty = ViennaFunctions::stack_energy(bp, next_bp, sequence) *
                               PseudoknotParams::pkp->pk_stack_x;
        if (round) stack_penalty = std::round(stack_penalty);
        return stack_penalty;
    }

    [[nodiscard]] static double pk_internal_energy(const BasePair& bp, const BasePair& next_bp,
                                                   const std::string& sequence, const bool& round) {
        double internal_penalty = ViennaFunctions::internal_loop_energy(bp, next_bp, sequence) *
                                  PseudoknotParams::pkp->pk_internal_x;
        if (round) internal_penalty = std::round(internal_penalty);
        return internal_penalty;
    }

    [[nodiscard]] static double pk_multiloop_energy(const BasePair& bp, const BasePair& next_bp,
                                                    const ProcessedRNAEntry& processed_rna) {
        double multiloop_penalty = PseudoknotParams::pkp->pk_mloop_init;

        // Since a multiloop is nested between two base pairs, we add 2 * bp_penalty
        // We add the child base pairs at a different part of the energy calculation
        multiloop_penalty += PseudoknotParams::pkp->pk_mloop_bp * 2;

        // get unpaired count
        int unpaired = processed_rna.get_unpaired_count(bp.i, next_bp.i);
        unpaired += processed_rna.get_unpaired_count(next_bp.j, bp.j);
        for (BasePair child_bp : bp.children) {
            unpaired -= processed_rna.get_unpaired_count(child_bp.i, child_bp.j);
        }

        // get unpaired penalty
        int pk_mloop_unpaired_energy = unpaired * PseudoknotParams::pkp->pk_mloop_unpaired;
        multiloop_penalty += pk_mloop_unpaired_energy;

        return multiloop_penalty;
    }
};

}  // namespace knotergy
