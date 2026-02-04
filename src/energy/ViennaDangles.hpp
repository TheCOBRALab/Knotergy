#pragma once

#include <array>
#include <vector>

#include "../io/ViennaParams.hpp"
#include "../loop_tree/LoopNode.hpp"
#include "ViennaUtils.hpp"

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

/**
 * @brief Stores dangle energy values for all four dangle configurations.
 *
 * Represents the energy contributions of dangling ends (unpaired nucleotides adjacent
 * to base pairs) in various configurations.
 */
struct DangleSet {
   public:
    DangleSet() : no_dangle(0), left_dangle(0), right_dangle(0), both_dangle(0) {}

    int no_dangle;      ///< Energy with no dangling ends.
    int left_dangle;    ///< Energy with only 5' dangling end.
    int right_dangle;   ///< Energy with only 3' dangling end.
    int both_dangle;    ///< Energy with both dangling ends.

    /**
     * @brief Get the minimum (most favorable) energy among all dangle configurations.
     *
     * @return The minimum energy value in centicalories.
     */
    int best() const { return std::min({no_dangle, left_dangle, right_dangle, both_dangle}); }

    /**
     * @brief Get the minimum energy between no dangle and left dangle only.
     *
     * @return The minimum energy value in centicalories.
     */
    int best_left() const { return std::min(no_dangle, left_dangle); }

    /**
     * @brief Get the minimum energy between no dangle and right dangle only.
     *
     * @return The minimum energy value in centicalories.
     */
    int best_right() const { return std::min(no_dangle, right_dangle); }
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
    static int get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children,
                                     const std::string& sequence);

    /**
     * @brief Calculate optimal external loop dangle energy using precomputed dangle sets.
     *
     * @param children Vector of child loop nodes in the external loop.
     * @param dangle_energies Precomputed dangle energies for each child.
     * @param sequence_length Length of the RNA sequence.
     * @return Optimal dangle energy in centicalories.
     */
    static int get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children,
                                     const std::vector<DangleSet>& dangle_energies, size_t sequence_length);

    /**
     * @brief Calculate optimal multibranch loop dangle energy.
     *
     * @param node The loop node representing the multibranch loop.
     * @param sequence The RNA nucleotide sequence.
     * @return Optimal dangle energy in centicalories.
     */
    static int get_multibranch_dangle_1(const LoopNode& node, const std::string& sequence);

    /**
     * @brief Calculate optimal multibranch loop dangle energy using precomputed values.
     *
     * @param node The loop node representing the multibranch loop.
     * @param dangle_energies Precomputed dangle energies for each child.
     * @param closing Dangle energy for the closing base pair.
     * @return Optimal dangle energy in centicalories.
     */
    static int get_multibranch_dangle_1(const LoopNode& node,
                                        std::vector<DangleSet> dangle_energies, DangleSet closing);

    /**
     * @brief Compute dangle energies for all child loop nodes.
     *
     * @param children Vector of child loop nodes.
     * @param sequence The RNA nucleotide sequence.
     * @param is_external Whether this is for an external loop (default: true).
     * @return Vector of DangleSet objects containing energies for each child.
     */
    static std::vector<DangleSet> populate_children_dangle_energies(
        const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence,
        const bool& is_external = true);

    /**
     * @brief Get dangle energy for a multibranch loop closing pair.
     *
     * @param node The loop node representing the multibranch loop.
     * @param sequence The RNA nucleotide sequence.
     * @return DangleSet containing energies for the closing pair.
     */
    static DangleSet get_ml_dangle_energy(const LoopNode& node, const std::string& sequence);

   private:
    /**
     * @brief Check if two indices are contiguous (adjacent).
     *
     * @param first First index.
     * @param second Second index.
     * @return True if indices differ by exactly 1.
     */
    static bool contiguous(size_t first, size_t second) noexcept {
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
    static bool contiguous_children(const LoopNode& first, const LoopNode& second) noexcept {
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
    static std::vector<std::vector<size_t>> get_dangle_chains(
        const std::vector<std::shared_ptr<LoopNode>>& children);

    /**
     * @brief Process a single chain to compute optimal dangle energy.
     *
     * @param chain Vector of child indices forming a contiguous chain.
     * @param children Vector of all child loop nodes.
     * @param dangle_energies Precomputed dangle energies.
     * @param disable_first_left_dangle Whether to disable left dangle for first child.
     * @param disable_last_right_dangle Whether to disable right dangle for last child.
     * @param init Initial values for dynamic programming.
     * @param closing Dangle energy for closing base pair (for multiloop).
     * @return Optimal dangle energy for the chain in centicalories.
     */
    static int process_chain(const std::vector<size_t>& chain,
                             const std::vector<std::shared_ptr<LoopNode>>& children,
                             const std::vector<DangleSet>& dangle_energies,
                             const bool& disable_first_left_dangle = false,
                             const bool& disable_last_right_dangle = false,
                             std::array<int, 2> init = {0, INF}, DangleSet closing = DangleSet());

    /**
     * @brief Process multiple chains to compute total optimal dangle energy.
     *
     * @param dangle_chains Vector of chains to process.
     * @param children Vector of all child loop nodes.
     * @param dangle_energies Precomputed dangle energies.
     * @param disable_first_left_dangle Whether to disable left dangle for first child.
     * @param disable_last_right_dangle Whether to disable right dangle for last child.
     * @param init Initial values for dynamic programming.
     * @param closing Dangle energy for closing base pair (for multiloop).
     * @return Total optimal dangle energy in centicalories.
     */
    static int process_chains(const std::vector<std::vector<size_t>>& dangle_chains,
                              const std::vector<std::shared_ptr<LoopNode>>& children,
                              const std::vector<DangleSet>& dangle_energies,
                              const bool& disable_first_left_dangle = false,
                              const bool& disable_last_right_dangle = false,
                              std::array<int, 2> init = {0, INF}, DangleSet closing = DangleSet());

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
    static int process_ml_chains(const std::vector<std::vector<size_t>>& dangle_chains,
                                 const std::vector<std::shared_ptr<LoopNode>>& children,
                                 const std::vector<DangleSet>& dangle_energies,
                                 const LoopNode& node, const DangleSet ml_dangle_energy);
};
}  // namespace knotergy