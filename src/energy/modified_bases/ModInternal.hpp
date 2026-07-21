#pragma once

#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/vienna/ViennaFunctions.hpp"

namespace knotergy {
class ModInternal {
   public:
    static int find_mod_internal_energy(size_t i, size_t j, size_t ci, size_t cj,
                                        const std::vector<std::string_view>& mod_sequence,
                                        const std::string& sequence, vrna_md_param& vp,
                                        const all_mod_params& mp);
};

}  // namespace knotergy