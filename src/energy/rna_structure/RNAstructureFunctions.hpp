#pragma once

#include "energy/rna_structure/RNAstructureUtils.hpp"
#include "io/parameters/RNAstructureParams.hpp"
#include "utils/common.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace knotergy {

/**
 * @brief efn2's (RNAstructure's) loop energies, computed from the tables in rnastructure_param.
 *
 * The functions follow RNAstructure's src/structure.cpp (erg1, erg3, ...) with 0-based
 * positions. Like ViennaFunctions, they return centicalories (hundredths of kcal/mol) and INF
 * for forbidden loops; RNAstructure's own tables are in tenths of kcal/mol.
 */
class RNAstructureFunctions {
   public:
    /**
     * @brief Calculate stacking energy for two consecutive base pairs (efn2's erg1).
     *
     * @param i 5' position of outer base pair.
     * @param j 3' position of outer base pair.
     * @param ci 5' position of inner base pair.
     * @param cj 3' position of inner base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views); empty without modified bases.
     * @param rsp RNAstructure parameters.
     * @return Stacking energy in centicalories.
     */
    [[nodiscard]] static int stack_energy(std::size_t i, std::size_t j, std::size_t ci,
                                          std::size_t cj, const std::string& sequence,
                                          const std::vector<std::string_view>& mod_sequence,
                                          const rnastructure_param& rsp);

    /**
     * @brief Calculate the helix-end term of a terminal base pair (efn2's penalty()).
     *
     * RNAstructure 6.6 replaced the terminal AU/GU penalty by the helix_ends table, which also
     * depends on the neighboring pair of the helix. The pair is given as seen from the loop:
     * the loop inside pair i-j (hairpin, internal, multibranch) uses helix_end(j, i), the loop
     * outside it uses helix_end(i, j). At the ends of the sequence the pair on the other side
     * is used, as in efn2.
     *
     * @param a Position of the pair that the loop meets first (5' to 3').
     * @param b Position paired with a.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views); empty without modified bases.
     * @param rsp RNAstructure parameters.
     * @return Helix-end energy in centicalories.
     */
    [[nodiscard]] static int helix_end(std::size_t a, std::size_t b, const std::string& sequence,
                                       const std::vector<std::string_view>& mod_sequence,
                                       const rnastructure_param& rsp);

    /**
     * @brief Calculate hairpin loop energy (efn2's erg3).
     *
     * Covers the loop initiation (extrapolated beyond 30 nucleotides), the first mismatch,
     * special tri-, tetra- and hexaloops, the GGG-U closure bonus, all-C loops, and the helix
     * end of the closing pair, which efn2 counts as part of the hairpin.
     *
     * @param i 5' position of closing base pair.
     * @param j 3' position of closing base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views); empty without modified bases.
     * @param rsp RNAstructure parameters.
     * @param is_inf Set to true when the energy is infinite (fewer than 3 loop nucleotides, a
     * closing pair RNAstructure does not allow, or a forbidden table entry such as N beside a
     * closing helix of one pair).
     * @return Hairpin loop energy in centicalories.
     */
    [[nodiscard]] static int hairpin_energy(std::size_t i, std::size_t j,
                                            const std::string& sequence,
                                            const std::vector<std::string_view>& mod_sequence,
                                            const rnastructure_param& rsp, bool& is_inf);

    /**
     * @brief Calculate internal loop or bulge energy (efn2's erg2).
     *
     * - A single bulged nucleotide: the two pairs stack across it, and the bulge gets a bonus
     *   of RT ln(n) for the n positions it can take among identical neighbors, plus the bonus
     *   of a C bulge next to another C. efn2.hpp adds the same two bonuses to ViennaRNA's
     *   energy; here they are already included.
     * - Longer bulges: initiation by size and the helix ends of both pairs.
     * - 1x1, 2x1 and 2x2 internal loops: their own tables.
     * - Other internal loops: initiation by size, the mismatches next to both pairs (with
     *   separate tables for 1xn and 2x3 loops) and the asymmetry penalty.
     * Internal loops also get the helix ends of both pairs. Sizes above 30 are extrapolated.
     * Without unpaired nucleotides the two pairs are a stack, and the stack energy is returned.
     *
     * @param i 5' position of outer base pair.
     * @param j 3' position of outer base pair.
     * @param ci 5' position of inner base pair.
     * @param cj 3' position of inner base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views); empty without modified bases.
     * @param rsp RNAstructure parameters.
     * @return Internal loop energy in centicalories; INF when RNAstructure does not allow the
     * loop: a pair it cannot form, or a forbidden table entry (N in a 1x1, 2x1 or 2x2 loop, or
     * beside a helix of one pair).
     */
    [[nodiscard]] static int internal_loop_energy(std::size_t i, std::size_t j, std::size_t ci,
                                                  std::size_t cj, const std::string& sequence,
                                                  const std::vector<std::string_view>& mod_sequence,
                                                  const rnastructure_param& rsp);
};

}  // namespace knotergy
