#pragma once

#include "io/parameters/ModParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
#include "utils/common.hpp"

namespace knotergy {

struct MultiloopStem {
    size_t begin = NULL_INDEX;  // Pairing base encountered while walking the multiloop
    size_t end = NULL_INDEX;    // Its pairing partner
    size_t prev_end =
        NULL_INDEX;  // The end of the previous stem in the walk, used for checking contiguity
    unsigned int type = 0;  // ViennaRNA pair type

    // Actual dangles
    int dangle5 = 0;  // Energy of 5' dangle for this stem
    int dangle3 = 0;  // Energy of 3' dangle for this stem

    // The initial 5' dangle energy (only first child and closing pair)
    int initial_ld5 = NULL_ENERGY;
};

class CoaxialStacking {
   public:
    /**
     * @brief Calculate multibranch loop dangle/coaxial contribution for dangle model 3.
     *
     * This uses a direct ViennaRNA-style fixed-structure d3 multiloop walk.
     * It is intentionally not implemented as dangle-1 plus a coaxial stacking bonus,
     * because ViennaRNA d3 uses separate state bookkeeping for odd dangles and
     * possible flush coaxial stacking.
     *
     * @param node The loop node representing the multibranch loop.
     * @param stems The multiloop stems in the loop, with precomputed dangle energies.
     * @param pRNA Processed RNA entry.
     * @param vp ViennaRNA model and parameter bundle.
     * @return Optimal dangle/coaxial contribution in centicalories.
     */
    [[nodiscard]] static int get_multibranch_dangle_3(const LoopNode& node,
                                                      const std::vector<MultiloopStem>& stems,
                                                      const ProcessedRNAEntry& pRNA,
                                                      vrna_md_param& vp,
                                                      const all_mod_params& mp = {});

    [[nodiscard]] static int get_multibranch_dangle_3(const LoopNode& node,
                                                      const ProcessedRNAEntry& pRNA,
                                                      vrna_md_param& vp,
                                                      const all_mod_params& mp = {});

    [[nodiscard]] static std::vector<MultiloopStem> populate_multiloop_stems(
        const LoopNode& node, const ProcessedRNAEntry& pRNA, vrna_md_param& vp);

   private:
    // ------------------- Dangle 3 (coaxial stacking) --------------------
    [[nodiscard]] static int compute_initial_ld5_for_d3(const MultiloopStem& stem,
                                                        const ProcessedRNAEntry& pRNA,
                                                        vrna_md_param& vp);
    [[nodiscard]] static int add_or_inf(int a, int b);

    [[nodiscard]] static int walk_multiloop_d3_from_start(const ProcessedRNAEntry& pRNA,
                                                          size_t start_prev,
                                                          const std::vector<MultiloopStem>& stems,
                                                          vrna_md_param& vp,
                                                          const all_mod_params& mp);
};

}  // namespace knotergy