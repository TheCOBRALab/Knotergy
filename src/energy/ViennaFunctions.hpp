#pragma once

#include "../loop_tree/LoopNode.hpp"
#include "../preprocessing/RNAEntry.hpp"
#include "ViennaDangles.hpp"
#include "ViennaUtils.hpp"
#include <algorithm>

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
class ViennaFunctions {
   public:
    ViennaFunctions(int dangle = 2) {
        vrna_md_set_default(&md);
        md.dangles = dangle;
        P = vrna_params(&md);
    }
    ~ViennaFunctions() { if (P) free(P); }

    const vrna_param_t* get_parameters() const;

    int stack_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence);
    int stack_energy(BasePair pair, BasePair child, const std::string& sequence);

    int hairpin_energy(size_t i, size_t j, const std::string& sequence);
    int hairpin_energy(const BasePair& pair, const std::string& sequence);

    int internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence);
    int internal_loop_energy(BasePair pair, BasePair child, const std::string& sequence);

    int multibranch_energy(const LoopNode& node, const std::string& sequence);

    int external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence);

   private:
    vrna_md_t md;
    vrna_param_t* P;
};
}  // namespace knotergy