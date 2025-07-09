#pragma once

#include <iostream>

#include "../loop_tree/LoopNode.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
class PseudoknotFunctions {
   public:
    PseudoknotFunctions() = default;
    ~PseudoknotFunctions() = default;
    // AGGGGUUUUUUUUUUUUUGGAAA
    // [[[[[.(((((]]]]]..)))))
    //
    double pseudoknot_energy(const LoopNode& node, const std::string& sequence) {
        const int P_tilda = 10;  // pair in a pseudoknot
        const int Q_tilda = 20;  // unpaired bases in a pseudoknot
        const int P_i = 10;      // for E&R energy model
        const int Gw = 700;      // starting a pseudoknot loop
        const int Gwh = 600;     // each more band region

        const double g_interiorPseudo = 0.83;
        [[maybe_unused]] const int multi_OffsetPseudo = 843;
        [[maybe_unused]] const int q_unpairedMultiPseudo = 0;
        [[maybe_unused]] const int p_pairedMultiPseudo = 100;

        double energy = 0;

        std::cout << node << std::endl;
        energy += Gw;
        energy += Gwh * (node.number_of_bands - 2);
        energy += P_tilda * 2 * node.number_of_bands;
        energy += Q_tilda * 2 * node.number_of_exclusive_unpaired_bases;
        energy += P_i * node.number_of_children_outside_band;

        for (const Band& band : node.bands) {
            const std::vector<Pair>& bp = band.base_pairs();
            for (size_t idx = 0; idx < bp.size() - 1; ++idx) {
                energy +=
                    g_interiorPseudo * vienna.stack_energy(bp[idx], bp[idx + 1], sequence) / 100.0;
            }
        }

        return energy;
    }

   private:
    ViennaFunctions vienna;
};
}  // namespace knotergy