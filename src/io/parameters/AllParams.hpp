#pragma once

#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/RNAstructureParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace knotergy {

enum class PKFreeEnergyModel {
    ViennaRNA2,  // ViennaRNA 2.x.x (Modified version of Turner 2004)
    efn2,        // RNAstructure's efn2 (Turner 2004)
    enf2Simple,  // RNAstructure's efn2 --simple (Turner 2004)
};

enum class PKEnergyModel {
    DP09,    // Dirks-Pierce 2009 (HotKnots V2)
    Hajdin,  // Hajdin 2013's polymer physics model (ShapeKnots/RNAstructure/efn2)
};

class AllParams {
   public:
    AllParams(const all_mod_params& mod, const pk_param& pk, const vrna_md_param& vienna,
              const rnastructure_param& rnastructure,
              PKFreeEnergyModel pk_free_energy_model = PKFreeEnergyModel::ViennaRNA2,
              PKEnergyModel pk_energy_model = PKEnergyModel::DP09)
        : mp_(mod),
          pkp_(pk),
          vp_(vienna),
          rsp_(rnastructure),
          pk_free_energy_model_(pk_free_energy_model),
          pk_energy_model_(pk_energy_model) {}

    const all_mod_params mp_;
    const pk_param pkp_;
    vrna_md_param vp_;
    const rnastructure_param rsp_;
    const PKFreeEnergyModel pk_free_energy_model_;
    const PKEnergyModel pk_energy_model_;
};

}  // namespace knotergy
