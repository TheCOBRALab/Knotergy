#pragma once
#include "io/parameters/AllParams.hpp"
#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/RNAstructureParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

namespace knotergy {

class EnergyFunctions {
   public:
    EnergyFunctions(AllParams& params)
        : params_(params),
          vp_(params.vp_),
          rsp_(params.rsp_),
          pkp_(params.pkp_),
          mp_(params.mp_),
          em_(params.pk_free_energy_model_),
          pk_em_(params.pk_energy_model_) {}

    [[nodiscard]] int stack(const std::size_t i, const std::size_t j, const std::size_t ci,
                            const std::size_t cj, const std::string& sequence,
                            const std::vector<std::string_view>& mod_sequence);
    [[nodiscard]] int d5(const std::size_t i, const std::size_t j, const std::string& sequence,
                         const std::vector<std::string_view>& mod_sequence,
                         const bool is_closing = false);
    [[nodiscard]] int d3(const std::size_t i, const std::size_t j, const std::string& sequence,
                         const std::vector<std::string_view>& mod_sequence,
                         const bool is_closing = false);
    [[nodiscard]] int mismatch(const std::size_t i, const std::size_t j,
                               const std::string& sequence,
                               const std::vector<std::string_view>& mod_sequence,
                               const bool is_external = false, const bool is_closing = false);
    [[nodiscard]] int terminal(const std::size_t i, const std::size_t j,
                               const std::string& sequence,
                               const std::vector<std::string_view>& mod_sequence,
                               const bool is_closing = false);
    [[nodiscard]] int hairpin(const std::size_t i, const std::size_t j, const std::string& sequence,
                              const std::vector<std::string_view>& mod_sequence, bool& is_inf);
    [[nodiscard]] int internal(const std::size_t i, const std::size_t j, const std::size_t ci,
                               const std::size_t cj, const std::string& sequence,
                               const std::vector<std::string_view>& mod_sequence);
    /// With RNAstructure's parameters: efn2's energy, or efn2 --simple's when
    /// simple_multibranch_rule is set.
    [[nodiscard]] int multibranch(const LoopNode& node, const ProcessedRNAEntry& pRNA);
    [[nodiscard]] int external(const std::vector<LoopNode*>& children,
                               const ProcessedRNAEntry& pRNA);

   private:
    const AllParams& params_;
    vrna_md_param& vp_;
    const rnastructure_param& rsp_;
    const pk_param& pkp_;
    const all_mod_params& mp_;
    const PKFreeEnergyModel em_;
    const PKEnergyModel pk_em_;
};
}  // namespace knotergy