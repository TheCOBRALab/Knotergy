#pragma once

#include "energy/dangles/DangleUtils.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

#include <algorithm>

namespace knotergy {

class Dangle1 {
   public:
    /**
     * @brief Calculate optimal external loop dangle energy for child loops.
     *
     * @param children Vector of child loop nodes in the external loop.
     * @param pRNA Processed RNA entry.
     * @param vp ViennaRNA model and parameter bundle.
     * @return Optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int get_external_dangle_1(const std::vector<LoopNode*>& children,
                                                   const ProcessedRNAEntry& pRNA,
                                                   vrna_md_param& vp);

    /**
     * @brief Calculate optimal external loop dangle energy using precomputed dangle sets.
     *
     * @param children Vector of child loop nodes in the external loop, in 5'->3' order.
     * @param dangle_energies Precomputed dangle energies for each child.
     * @return Optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int get_external_dangle_1(const std::vector<LoopNode*>& children,
                                                   const std::vector<DangleSet>& dangle_energies);

    /**
     * @brief Calculate optimal multibranch loop dangle energy.
     *
     * @param node The loop node representing the multibranch loop.
     * @param pRNA Processed RNA entry.
     * @param vp ViennaRNA model and parameter bundle.
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
     * @param closing Dangle energies of the closing pair, seen from inside the loop.
     * @return Optimal dangle energy in centicalories.
     */
    [[nodiscard]] static int get_multibranch_dangle_1(const LoopNode& node,
                                                      const std::vector<DangleSet>& dangle_energies,
                                                      const DangleSet& closing);

    /**
     * @brief Compute dangle energies for all child loop nodes.
     *
     * @param children Vector of child loop nodes.
     * @param pRNA Processed RNA entry.
     * @param vp ViennaRNA model and parameter bundle.
     * @param is_external Whether this is for an external loop (default: true).
     * @return Vector of DangleSet objects containing energies for each child.
     */
    [[nodiscard]] static std::vector<DangleSet> populate_children_dangle_energies(
        const std::vector<LoopNode*>& children, const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
        bool is_external = true);

    /**
     * @brief Get dangle energy for a multibranch loop closing pair.
     *
     * @param node The loop node representing the multibranch loop.
     * @param pRNA Processed RNA entry.
     * @param vp ViennaRNA model and parameter bundle.
     * @return DangleSet containing energies for the closing pair.
     */
    [[nodiscard]] static DangleSet get_ml_closing_dangle_energy(const LoopNode& node,
                                                                const ProcessedRNAEntry& pRNA,
                                                                vrna_md_param& vp);

    [[nodiscard]] static DangleSet get_child_dangle_energy(const LoopNode& node,
                                                           const ProcessedRNAEntry& pRNA,
                                                           vrna_md_param& vp,
                                                           bool is_external = true);
};
};  // namespace knotergy
