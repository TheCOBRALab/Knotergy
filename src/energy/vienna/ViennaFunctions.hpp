#pragma once

#include "energy/dangles/CoaxialStacking.hpp"
#include "energy/dangles/Dangle1.hpp"
#include "energy/vienna/ViennaUtils.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/RNAEntry.hpp"

#include <algorithm>

namespace knotergy {
/**
 * @brief Provides ViennaRNA-based energy calculation functions for RNA secondary structures.
 *
 * This class wraps ViennaRNA library functions to compute energies for different loop types:
 * stacking pairs, hairpin loops, internal loops, multiloop branches, and external loops.
 * All energy values are returned in centicalories (hundredths of kcal/mol).
 */
class ViennaFunctions {
   public:
    /**
     * @brief Calculate stacking energy for two consecutive base pairs.
     *
     * @param i 5' position of outer base pair.
     * @param j 3' position of outer base pair.
     * @param ci 5' position of inner base pair.
     * @param cj 3' position of inner base pair.
     * @param sequence The RNA nucleotide sequence.
     * @return Stacking energy in centicalories.
     */
    [[nodiscard]] static int stack_energy(size_t i, size_t j, size_t ci, size_t cj,
                                          const std::string& sequence, vrna_md_param& vp);

    /**
     * @brief Calculate stacking energy for two consecutive base pairs.
     *
     * @param node The loop node representing the stacking pair.
     * @param sequence The RNA nucleotide sequence.
     * @param vp The ViennaRNA model parameters.
     * @return Stacking energy in centicalories.
     */
    [[nodiscard]] static int stack_energy(const LoopNode& node, const std::string& sequence,
                                          vrna_md_param& vp);

    /**
     * @brief Calculate stacking energy for two consecutive base pairs.
     *
     * @param pair The outer base pair.
     * @param child The inner base pair.
     * @param sequence The RNA nucleotide sequence.
     * @return Stacking energy in centicalories.
     */
    [[nodiscard]] static int stack_energy(BasePair pair, BasePair child,
                                          const std::string& sequence, vrna_md_param& vp);

    /**
     * @brief Calculate hairpin loop energy.
     *
     * @param i 5' position of closing base pair.
     * @param j 3' position of closing base pair.
     * @param sequence The RNA nucleotide sequence.
     * @param is_inf Whether the energy is infinite (hairpin loop size < 3).
     * @return Hairpin loop energy in centicalories.
     */
    [[nodiscard]] static int hairpin_energy(size_t i, size_t j, const std::string& sequence,
                                            bool& is_inf, vrna_md_param& vp);

    /**
     * @brief Calculate hairpin loop energy.
     *
     * @param node The loop node representing the hairpin loop.
     * @param sequence The RNA nucleotide sequence.
     * @param is_inf Whether the energy is infinite (hairpin loop size < 3).
     * @return Hairpin loop energy in centicalories.
     */
    [[nodiscard]] static int hairpin_energy(const LoopNode& node, const std::string& sequence,
                                            bool& is_inf, vrna_md_param& vp);

    /**
     * @brief Calculate hairpin loop energy.
     *
     * @param pair The closing base pair of the hairpin loop.
     * @param sequence The RNA nucleotide sequence.
     * @param is_inf Whether the energy is infinite (hairpin loop size < 3).
     * @return Hairpin loop energy in centicalories.
     */
    [[nodiscard]] static int hairpin_energy(const BasePair& pair, const std::string& sequence,
                                            bool& is_inf, vrna_md_param& vp);

    /**
     * @brief Calculate internal loop or bulge energy.
     *
     * @param i 5' position of outer base pair.
     * @param j 3' position of outer base pair.
     * @param ci 5' position of inner base pair.
     * @param cj 3' position of inner base pair.
     * @param sequence The RNA nucleotide sequence.
     * @return Internal loop energy in centicalories.
     */
    [[nodiscard]] static int internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj,
                                                  const std::string& sequence, vrna_md_param& vp);

    [[nodiscard]] static int internal_loop_energy(const LoopNode& node, const std::string& sequence,
                                                  vrna_md_param& vp);

    /**
     * @brief Calculate internal loop or bulge energy.
     *
     * @param pair The outer base pair.
     * @param child The inner base pair.
     * @param sequence The RNA nucleotide sequence.
     * @return Internal loop energy in centicalories.
     */
    [[nodiscard]] static int internal_loop_energy(BasePair pair, BasePair child,
                                                  const std::string& sequence, vrna_md_param& vp);

    /**
     * @brief Calculate multibranch loop energy.
     *
     * @param node The loop node representing the multibranch loop.
     * @param sequence The RNA nucleotide sequence.
     * @return Multibranch loop energy in centicalories.
     */
    [[nodiscard]] static int multibranch_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                                vrna_md_param& vp);

    /**
     * @brief Calculate external loop energy.
     *
     * @param children Vector of child loop nodes in the external loop.
     * @param sequence The RNA nucleotide sequence.
     * @return External loop energy in centicalories.
     */
    [[nodiscard]] static int external_energy(const std::vector<std::unique_ptr<LoopNode>>& children,
                                             const ProcessedRNAEntry& pRNA, vrna_md_param& vp);
};
}  // namespace knotergy