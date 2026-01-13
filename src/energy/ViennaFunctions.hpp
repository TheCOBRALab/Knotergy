#pragma once

#include "../loop_tree/LoopNode.hpp"
#include "../preprocessing/RNAEntry.hpp"
#include "../pipeline/load_params.hpp"
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
    static int stack_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence);
    static int stack_energy(BasePair pair, BasePair child, const std::string& sequence);

    static int hairpin_energy(size_t i, size_t j, const std::string& sequence);
    static int hairpin_energy(const BasePair& pair, const std::string& sequence);

    static int internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence);
    static int internal_loop_energy(BasePair pair, BasePair child, const std::string& sequence);

    static int multibranch_energy(const LoopNode& node, const std::string& sequence);

    static int external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence);

};
}  // namespace knotergy