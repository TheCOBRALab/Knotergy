#include "io/parameters/AllParams.hpp"

#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/modified_bases/ModHairpin.hpp"
#include "energy/modified_bases/ModInternal.hpp"
#include "energy/modified_bases/ModStack.hpp"
#include "energy/rna_structure/RNAstructureFunctions.hpp"
#include "energy/vienna/ViennaFunctions.hpp"
#include "utils/common.hpp"

namespace knotergy {
int AllParams::stack(const std::size_t i, const std::size_t j, const std::size_t ci,
                     const std::size_t cj, const std::string& sequence,
                     const std::vector<std::string_view>& mod_sequence) {
    if (rnastructure_params_loaded) {
        return RNAstructureFunctions::stack_energy(i, j, ci, cj, sequence, mod_sequence, rsp_);
    }

    if (modified_params_loaded) {
        return ModStack::find_mod_stack_energy(i, j, ci, cj, sequence, mod_sequence, vp_, mp_);
    }

    if (vienna_params_loaded) {
        return ViennaFunctions::stack_energy(i, j, ci, cj, sequence, vp_);
    }

    return NULL_ENERGY;
}

int AllParams::d5(const std::size_t i, const std::size_t j, const std::string& sequence,
                  const std::vector<std::string_view>& mod_sequence, const bool is_closing) {
    if (rnastructure_params_loaded) {
        return RNAstructureFunctions::dangle5_energy(i, j, sequence, mod_sequence, rsp_,
                                                     is_closing);
    }

    if (modified_params_loaded) {
        // The JSON files only have the dangles that involve a modified base.
        const int energy = ModBaseUtils::danlge5_mod_energy_lookup(i, j, mod_sequence, mp_,
                                                                   vp_.md.temperature, is_closing);
        return energy != NULL_ENERGY
                   ? energy
                   : ViennaFunctions::dangle5_energy(i, j, sequence, vp_, is_closing);
    }

    if (vienna_params_loaded) {
        return ViennaFunctions::dangle5_energy(i, j, sequence, vp_, is_closing);
    }

    return NULL_ENERGY;
}

int AllParams::d3(const std::size_t i, const std::size_t j, const std::string& sequence,
                  const std::vector<std::string_view>& mod_sequence, const bool is_closing) {
    if (rnastructure_params_loaded) {
        return RNAstructureFunctions::dangle3_energy(i, j, sequence, mod_sequence, rsp_,
                                                     is_closing);
    }

    if (modified_params_loaded) {
        // The JSON files only have the dangles that involve a modified base.
        const int energy = ModBaseUtils::danlge3_mod_energy_lookup(i, j, mod_sequence, mp_,
                                                                   vp_.md.temperature, is_closing);
        return energy != NULL_ENERGY
                   ? energy
                   : ViennaFunctions::dangle3_energy(i, j, sequence, vp_, is_closing);
    }

    if (vienna_params_loaded) {
        return ViennaFunctions::dangle3_energy(i, j, sequence, vp_, is_closing);
    }

    return NULL_ENERGY;
}

int AllParams::hairpin(const std::size_t i, const std::size_t j, const std::string& sequence,
                       const std::vector<std::string_view>& mod_sequence, bool& is_inf) {
    if (rnastructure_params_loaded) {
        return RNAstructureFunctions::hairpin_energy(i, j, sequence, mod_sequence, rsp_, is_inf);
    }

    if (modified_params_loaded) {
        return ModHairpin::find_mod_hairpin_energy(i, j, sequence, mod_sequence, vp_, mp_, is_inf);
    }

    if (vienna_params_loaded) {
        return ViennaFunctions::hairpin_energy(i, j, sequence, is_inf, vp_);
    }

    return NULL_ENERGY;
}

int AllParams::internal(const std::size_t i, const std::size_t j, const std::size_t ci,
                        const std::size_t cj, const std::string& sequence,
                        const std::vector<std::string_view>& mod_sequence) {
    if (rnastructure_params_loaded) {
        return RNAstructureFunctions::internal_loop_energy(i, j, ci, cj, sequence, mod_sequence,
                                                           rsp_);
    }

    if (modified_params_loaded) {
        return ModInternal::find_mod_internal_energy(i, j, ci, cj, mod_sequence, sequence, vp_,
                                                     mp_);
    }

    if (vienna_params_loaded) {
        return ViennaFunctions::internal_loop_energy(i, j, ci, cj, sequence, vp_);
    }

    return NULL_ENERGY;
}

}  // namespace knotergy
