#pragma once

#include <iostream>

#include "../loops/LoopNode.hpp"

namespace knotergy {
class PseudoknotFunctions {
   public:
    PseudoknotFunctions() = default;
    ~PseudoknotFunctions() = default;
    // AGGGGUUUUUUUUUUUUUGGAAA
    // [[[[[.(((((]]]]]..)))))
    //
    float pseudoknot_energy(const LoopNode& node) {
        double P_tilda = 0.1f;  // pair in a pseudoknot
        double Q_tilda = 0.2f;  // unpaired bases in a pseudoknot
        double P_i = 0.1f;      // for E&R energy model
        double Gw = 7.0f;       // starting a pseudoknot loop
        double Gwh = 6.0f;      // each more band region

        double energy = 0.0f;

        std::cout << node << std::endl;
        energy += Gw;
        energy += Gwh * (node.number_of_bands - 2);
        energy += P_tilda * 2 * node.number_of_bands;
        energy += Q_tilda * 2 * node.number_of_exclusive_unpaired_bases;
        energy += P_i * node.number_of_children_outside_band;
        energy *= 100;

        return static_cast<float>(energy);
    }

   private:
};
}  // namespace knotergy