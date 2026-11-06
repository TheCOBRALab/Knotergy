#include "io/parameters/AllParams.hpp"

#include "energy/modified_bases/ModStack.hpp"
#include "energy/rna_structure/RNAstructureFunctions.hpp"
#include "energy/vienna/ViennaFunctions.hpp"
#include "utils/common.hpp"

namespace knotergy {
int AllParams::stack(const std::size_t i, const std::size_t j, const std::size_t ci,
                     const std::size_t cj, const std::string& sequence,
                     const std::vector<std::string_view>& mod_sequence) {
    int energy = NULL_ENERGY;

    if (modified_params_loaded) {
        energy = ModStack::find_mod_stack_energy(i, j, ci, cj, sequence, mod_sequence, vp_, mp_);
    }

    if (energy == NULL_ENERGY && rnastructure_params_loaded) {
        return RNAstructureFunctions::stack_energy(i, j, ci, cj, sequence, mod_sequence, rsp_);
    }

    if (vienna_params_loaded) {
        return ViennaFunctions::stack_energy(i, j, ci, cj, sequence, vp_);
    }

    return energy;
}

}  // namespace knotergy