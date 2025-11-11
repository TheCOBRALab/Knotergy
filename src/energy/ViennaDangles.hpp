#pragma once

#include <vector>
#include <array>
#include "../loop_tree/LoopNode.hpp"
#include "ViennaUtils.hpp"

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

    struct DangleSet {
    public:
        DangleSet() : no_dangle(0), left_dangle(0), right_dangle(0), both_dangle(0) {}
        int no_dangle;
        int left_dangle;
        int right_dangle;
        int both_dangle;

        int min() const {
            return std::min({no_dangle, left_dangle, right_dangle, both_dangle});
        }
    };

    class ViennaDangles {
    public:
        ViennaDangles() = default;
        ~ViennaDangles() = default;
        static int get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence, vrna_md_t& md);
        static int get_multi_dangle_1(const LoopNode& node, const std::string& sequence, vrna_md_t& md);
    
    private:
        
        static std::vector<DangleSet> populate_children_dangle_energies(
            const std::vector<std::shared_ptr<LoopNode>>& children,
            const std::string& sequence,
            vrna_md_t& md,
            const bool& is_external = true
        );
        static std::vector<std::vector<size_t>> get_dangle_chains(const std::vector<std::shared_ptr<LoopNode>>& children);
        static int process_chain(
                    const std::vector<size_t>& chain,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const bool& disable_first_left_dangle = false,
                    const bool& disable_first_right_dangle = false
                );
        static int process_chains(
                    const std::vector<std::vector<size_t>>& dangle_chains,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies
        );

        static DangleSet get_ml_dangle_energy(const LoopNode& node, const std::string& sequence, vrna_md_t& md);
        static int process_ml_chains(
                    const std::vector<std::vector<size_t>>& dangle_chains,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const LoopNode& node,
                    const DangleSet ml_dangle_energy
                );
    };
} // namespace knotergy