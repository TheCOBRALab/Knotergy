#pragma once

#include "energy/rna_structure/RNAstructureUtils.hpp"
#include "io/parameters/RNAstructureParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
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

    /**
     * @brief Calculate the energy of a 5' dangling end on base pair i-j (efn2's erg4).
     *
     * A dangling end is an unpaired nucleotide of an exterior or multibranch loop that stacks
     * on the end of a helix. Going through the loop 5' to 3', the 5' dangle is the nucleotide
     * just before the helix: i - 1 for a branch of the loop, j - 1 for the pair that closes a
     * multibranch loop (is_closing), whose dangles are inside it.
     *
     * efn2 uses a dangle only where it gives the lowest energy among the ways the nucleotides
     * next to a helix can stack (5' dangle, 3' dangle, terminal mismatch, coaxial stack, or
     * none of them), so this is a term to choose from, not one to always add.
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views); empty without modified bases.
     * @param rsp RNAstructure parameters.
     * @param is_closing True when i-j closes a multibranch loop: the dangle is inside the pair.
     * @return Dangle energy in centicalories; INF where RNAstructure's table has no value: a
     * pair it cannot form, or N, which it does not allow as a 5' dangle.
     */
    [[nodiscard]] static int dangle5_energy(std::size_t i, std::size_t j,
                                            const std::string& sequence,
                                            const std::vector<std::string_view>& mod_sequence,
                                            const rnastructure_param& rsp, bool is_closing = false);

    /**
     * @brief Calculate the energy of a 3' dangling end on base pair i-j (efn2's erg4).
     *
     * Going through the loop 5' to 3', the 3' dangle is the nucleotide just after the helix:
     * j + 1 for a branch of an exterior or multibranch loop, i + 1 for the pair that closes a
     * multibranch loop (is_closing). See dangle5_energy.
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views); empty without modified bases.
     * @param rsp RNAstructure parameters.
     * @param is_closing True when i-j closes a multibranch loop: the dangle is inside the pair.
     * @return Dangle energy in centicalories, 0 when the dangling nucleotide is N; otherwise INF
     * for a pair RNAstructure cannot form.
     */
    [[nodiscard]] static int dangle3_energy(std::size_t i, std::size_t j,
                                            const std::string& sequence,
                                            const std::vector<std::string_view>& mod_sequence,
                                            const rnastructure_param& rsp, bool is_closing = false);

    /**
     * @brief Calculate the terminal mismatch energy of base pair i-j in an exterior or
     * multibranch loop.
     *
     * A terminal mismatch is the two unpaired nucleotides on either side of the end of a
     * helix, stacked on it together: i - 1 and j + 1 for a branch of the loop, j - 1 and i + 1
     * for the pair that closes a multibranch loop (is_closing). efn2 has one table for exterior
     * loops (is_external) and one for multibranch loops. Like a dangle, it is a term to choose
     * from (see dangle5_energy), and it does not include the helix end (see helix_end).
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views); empty without modified bases.
     * @param rsp RNAstructure parameters.
     * @param is_external True in an exterior loop, false in a multibranch loop.
     * @param is_closing True when i-j closes a multibranch loop: the mismatch is inside the pair.
     * @return Mismatch energy in centicalories; INF for a pair RNAstructure cannot form.
     */
    [[nodiscard]] static int mismatch_energy(std::size_t i, std::size_t j,
                                             const std::string& sequence,
                                             const std::vector<std::string_view>& mod_sequence,
                                             const rnastructure_param& rsp,
                                             bool is_external = false, bool is_closing = false);

    /**
     * @brief Calculate multibranch loop energy (efn2's ergmulti).
     *
     * The loop is a circle of helices with unpaired nucleotides between them: the pair that
     * closes it, seen from inside, and its branches. Its energy is the sum of
     * - an offset, a penalty per helix and a penalty per unpaired nucleotide, which is
     *   logarithmic beyond 8 nucleotides;
     * - an asymmetry penalty, for how much the numbers of unpaired nucleotides on the two sides
     *   of a helix differ, averaged over the helices;
     * - a strain penalty when an odd number of helices has at most one unpaired nucleotide
     *   among them (a three-way junction without unpaired nucleotides, for example);
     * - the helix end of every helix;
     * - the lowest energy over the ways the helices and the nucleotides next to them can stack,
     *   as in an exterior loop (see external_energy), with the mismatches of multibranch loops.
     *
     * This is efn2's default. Its --simple option, which leaves out the logarithm, the
     * asymmetry and the strain like RNAstructure's folding algorithms do, is not implemented.
     *
     * A pseudoknotted child counts as a helix of the loop, but it has no helix end and nothing
     * stacks on it, as in an exterior loop.
     *
     * @param node The loop node representing the multibranch loop.
     * @param pRNA Processed RNA entry.
     * @param rsp RNAstructure parameters.
     * @return Multibranch loop energy in centicalories; INF when RNAstructure does not allow the
     * end of one of the helices (a pair it cannot form, for example).
     */
    [[nodiscard]] static int multibranch_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                                const rnastructure_param& rsp);

    /**
     * @brief Calculate exterior loop energy (efn2's ergexterior).
     *
     * The helix end of every branch, plus the lowest energy over the ways the helices and the
     * nucleotides next to them can stack. Going through the loop 5' to 3', a helix does at most
     * one of these, and a nucleotide is used at most once:
     * - take the nucleotide before it as a 5' dangle, or the one after it as a 3' dangle;
     * - take both as a terminal mismatch;
     * - stack coaxially on the next helix when no nucleotide separates them;
     * - stack coaxially on the next helix across the single nucleotide between them, which
     *   forms a mismatch with the nucleotide before the first helix or after the second.
     *
     * A pseudoknotted child is not a helix RNAstructure has stacking rules for: nothing stacks
     * on it, and the branches on its two sides do not stack on each other.
     *
     * @param children Child loop nodes of the exterior loop, in 5' to 3' order.
     * @param pRNA Processed RNA entry.
     * @param rsp RNAstructure parameters.
     * @return Exterior loop energy in centicalories; INF when RNAstructure does not allow the
     * end of a branch (a pair it cannot form, for example).
     */
    [[nodiscard]] static int external_energy(const std::vector<LoopNode*>& children,
                                             const ProcessedRNAEntry& pRNA,
                                             const rnastructure_param& rsp);
};

}  // namespace knotergy
