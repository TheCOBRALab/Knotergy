#pragma once

#include "energy/vienna/ViennaUtils.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

extern "C" {
#include <ViennaRNA/eval/exterior.h>
#include <ViennaRNA/eval/hairpin.h>
#include <ViennaRNA/eval/internal.h>
#include <ViennaRNA/eval/multibranch.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/sequences/alphabet.h>
#include <ViennaRNA/utils/basic.h>
}

namespace knotergy {

struct MultiloopStem {
    size_t begin;       // Pairing base encountered while walking the multiloop
    size_t end;         // Its pairing partner
    size_t prev_end;    // The end of the previous stem in the walk, used for checking contiguity
    unsigned int type;  // ViennaRNA pair type

    // Actual dangles
    int dangle5;  // Energy of 5' dangle for this stem
    int dangle3;  // Energy of 3' dangle for this stem

    // The initial 5' dangle energy (only first child and closing pair)
    int initial_ld5 = NULL_ENERGY;
};

// ------------------ Dangle 1 --------------------

/**
 * @brief Stores dangle energy values for all four dangle configurations.
 *
 * Represents the energy contributions of dangling ends (unpaired nucleotides adjacent
 * to base pairs) in various configurations.
 */
struct DangleSet {
   public:
    DangleSet() : no_dangle(0), left_dangle(0), right_dangle(0), both_dangle(0) {}
    DangleSet(int no_dangle_energy, int left_dangle_energy, int right_dangle_energy,
              int both_dangle_energy)
        : no_dangle(no_dangle_energy),
          left_dangle(left_dangle_energy),
          right_dangle(right_dangle_energy),
          both_dangle(both_dangle_energy) {}

    int no_dangle;     ///< Energy with no dangling ends.
    int left_dangle;   ///< Energy with only 5' dangling end.
    int right_dangle;  ///< Energy with only 3' dangling end.
    int both_dangle;   ///< Energy with both dangling ends.

    /**
     * @brief Get the minimum (most favorable) energy among all dangle configurations.
     *
     * @return The minimum energy value in centicalories.
     */
    [[nodiscard]] int best() const {
        return std::min({no_dangle, left_dangle, right_dangle, both_dangle});
    }

    /**
     * @brief Get the minimum energy between no dangle and left dangle only.
     *
     * @return The minimum energy value in centicalories.
     */
    [[nodiscard]] int best_left() const { return std::min(no_dangle, left_dangle); }

    /**
     * @brief Get the minimum energy between no dangle and right dangle only.
     *
     * @return The minimum energy value in centicalories.
     */
    [[nodiscard]] int best_right() const { return std::min(no_dangle, right_dangle); }

    /** @brief Add a constant energy to all dangle configurations.
     *
     * This is useful for applying energy contributions that affect all configurations
     * equally, such as terminal mismatches
     *
     * @param energy The energy value in centicalories to add to each configuration.
     * @return Reference to the modified DangleSet.
     */
    DangleSet& operator+=(int energy) {
        no_dangle += energy;
        left_dangle += energy;
        right_dangle += energy;
        both_dangle += energy;
        return *this;
    }
};

/**
 * @brief Handles dangling end energy calculations for RNA secondary structures.
 *
 * Implements the dangle model for calculating the stabilizing effect of unpaired
 * nucleotides adjacent to base pairs. Uses dynamic programming to find optimal
 * dangle configurations for external and multibranch loops.
 */
class ViennaDangles {
   public:
    ViennaDangles() = default;
    ~ViennaDangles() = default;

    /**
     * @brief Calculate optimal external loop dangle energy for child loops.
     *
     * @param children Vector of child loop nodes in the external loop.
     * @param sequence The RNA nucleotide sequence.
     * @return Optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int get_external_dangle_1(
        const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
        vrna_md_param& vp);

    /**
     * @brief Calculate optimal external loop dangle energy using precomputed dangle sets.
     *
     * @param children Vector of child loop nodes in the external loop.
     * @param dangle_energies Precomputed dangle energies for each child.
     * @param sequence_length Length of the RNA sequence.
     * @return Optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int get_external_dangle_1(
        const std::vector<std::unique_ptr<LoopNode>>& children,
        const std::vector<DangleSet>& dangle_energies);

    /**
     * @brief Calculate optimal multibranch loop dangle energy.
     *
     * @param node The loop node representing the multibranch loop.
     * @param sequence The RNA nucleotide sequence.
     * @return Optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int get_multibranch_dangle_1(const LoopNode& node,
                                                      const ProcessedRNAEntry& pRNA,
                                                      vrna_md_param& vp);

    /**
     * @brief Calculate optimal multibranch loop dangle energy using precomputed values.
     *
     * @param node The loop node representing the multibranch loop.
     * @param dangle_energies Precomputed dangle energies for each child.
     * @param closing Dangle energy for the closing base pair.
     * @return Optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int get_multibranch_dangle_1(const LoopNode& node,
                                                      std::vector<DangleSet> dangle_energies,
                                                      DangleSet closing);

    /**
     * @brief Calculate multibranch loop dangle/coaxial contribution for dangle model 3.
     *
     * This uses a direct ViennaRNA-style fixed-structure d3 multiloop walk.
     * It is intentionally not implemented as dangle-1 plus a coaxial stacking bonus,
     * because ViennaRNA d3 uses separate state bookkeeping for odd dangles and
     * possible flush coaxial stacking.
     *
     * @param node The loop node representing the multibranch loop.
     * @param pRNA Processed RNA entry.
     * @param vp ViennaRNA model and parameter bundle.
     * @return Optimal dangle/coaxial contribution in centicalories.
     */
    [[nodiscard]] static int get_multibranch_dangle_3(const LoopNode& node,
                                                      const ProcessedRNAEntry& pRNA,
                                                      vrna_md_param& vp);

    /**
     * @brief Compute dangle energies for all child loop nodes.
     *
     * @param children Vector of child loop nodes.
     * @param sequence The RNA nucleotide sequence.
     * @param is_external Whether this is for an external loop (default: true).
     * @return Vector of DangleSet objects containing energies for each child.
     */
    [[nodiscard]] static std::vector<DangleSet> populate_children_dangle_energies(
        const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
        vrna_md_param& vp, bool is_external = true);

    /**
     * @brief Get dangle energy for a multibranch loop closing pair.
     *
     * @param node The loop node representing the multibranch loop.
     * @param sequence The RNA nucleotide sequence.
     * @return DangleSet containing energies for the closing pair.
     */
    [[nodiscard]] static DangleSet get_ml_closing_dangle_energy(const LoopNode& node,
                                                                const ProcessedRNAEntry& pRNA,
                                                                vrna_md_param& vp);

    [[nodiscard]] static DangleSet get_child_dangle_energy(const LoopNode& node,
                                                           const ProcessedRNAEntry& pRNA,
                                                           vrna_md_param& vp,
                                                           bool is_external = true);

   private:
    /**
     * @brief Check if two indices are contiguous (adjacent).
     *
     * @param first First index.
     * @param second Second index.
     * @return True if indices differ by exactly 1.
     */
    [[nodiscard]] static bool contiguous(size_t first, size_t second) noexcept {
        return (first > second ? first - second : second - first) == 1;
    }

    /**
     * @brief Check if two child loop nodes are directly adjacent in sequence.
     *
     * Assumes children are in 5'→3' order but handles reversed input defensively.
     *
     * @param first First loop node.
     * @param second Second loop node.
     * @return True if the nodes are adjacent in the sequence.
     */
    [[nodiscard]] static bool contiguous_children(const LoopNode& first,
                                                  const LoopNode& second) noexcept {
        if (first.end < second.begin) return contiguous(first.end, second.begin);
        if (second.end < first.begin) return contiguous(second.end, first.begin);
        return false;  // overlapping or nested, not adjacent
    }

    /**
     * @brief Identify chains of contiguous children for dangle processing.
     *
     * @param children Vector of child loop nodes.
     * @return Vector of chains, where each chain is a vector of child indices.
     */
    [[nodiscard]] static std::vector<std::vector<size_t>> get_dangle_chains(
        const std::vector<std::unique_ptr<LoopNode>>& children);

    /**
     * @brief Process a single chain to compute optimal dangle energy.
     *
     * @param chain Vector of child indices forming a contiguous chain.
     * @param children Vector of all child loop nodes.
     * @param dangle_energies Precomputed dangle energies.
     * @param disable_last_right_dangle Whether to disable right dangle for last child.
     * @param init Initial values for dynamic programming.
     * @param closing Dangle energy for closing base pair (for multiloop).
     * @return Optimal dangle energy for the chain in centicalories.
     */
    [[nodiscard]] static int process_chain(const std::vector<size_t>& chain,
                                           const std::vector<DangleSet>& dangle_energies,
                                           bool disable_last_right_dangle = false,
                                           std::array<int, 2> init = {0, INF},
                                           DangleSet closing = DangleSet());

    /**
     * @brief Process multiple chains to compute total optimal dangle energy.
     *
     * @param dangle_chains Vector of chains to process.
     * @param children Vector of all child loop nodes.
     * @param dangle_energies Precomputed dangle energies.
     * @param disable_last_right_dangle Whether to disable right dangle for last child.
     * @param init Initial values for dynamic programming.
     * @param closing Dangle energy for closing base pair (for multiloop).
     * @return Total optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int process_chains(const std::vector<std::vector<size_t>>& dangle_chains,
                                            const std::vector<DangleSet>& dangle_energies,
                                            bool disable_last_right_dangle = false,
                                            std::array<int, 2> init = {0, INF},
                                            DangleSet closing = DangleSet());

    /**
     * @brief Process chains for multibranch loop dangle energy calculation.
     *
     * @param dangle_chains Vector of chains to process.
     * @param children Vector of all child loop nodes.
     * @param dangle_energies Precomputed dangle energies.
     * @param node The multibranch loop node.
     * @param ml_dangle_energy Dangle energy for the multibranch loop closing pair.
     * @return Total optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int process_ml_chains(
        const std::vector<std::vector<size_t>>& dangle_chains,
        const std::vector<std::unique_ptr<LoopNode>>& children,
        const std::vector<DangleSet>& dangle_energies, const LoopNode& node,
        const DangleSet ml_dangle_energy);

    // --------------- Dangle 3: Coaxial Stacking ------------------------
};
}  // namespace knotergy