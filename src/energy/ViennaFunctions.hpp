#pragma once

#include "../loop_tree/LoopNode.hpp"
#include "../preprocessing/RNAEntry.hpp"

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
        std::cout << "Initialized ViennaFunctions with dangle model: " << dangle << std::endl;
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
    enum DangleIdx { None = 0, Left = 1, Right = 2, Both = 3 };

    unsigned int get_pair_type(const char& i, const char& j);

    unsigned int reverse_pair_type(unsigned int type) const;

    int get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence);
    int get_multi_dangle_1(const LoopNode& node, const std::string& sequence);
    std::vector<std::array<int,4>> populate_children_dangle_energies(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence, const bool& is_external=true);
    std::array<int,4> populate_ml_dangle_energies(const LoopNode& node, const std::string& sequence);
    std::vector<std::vector<size_t>> get_dangle_chains(const std::vector<std::shared_ptr<LoopNode>>& children);
};
}  // namespace knotergy