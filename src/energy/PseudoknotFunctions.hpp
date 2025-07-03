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
    float pseudoknot_energy(const LoopNode& node) {
        float P_tilda = 0.1f;  // pair in a pseudoknot
        float Q_tilda = 0.2f;  // unpaired bases in a pseudoknot
        float P_i = 0.1f;      // for E&R energy model
        float Gw = 7.0f;       // starting a pseudoknot loop
        float Gwh = 6.0f;      // each more band region

        float energy = 0.0f;

        std::cout << "Band Count: " << node.number_of_bands
                  << " Unpaired Count: " << node.number_of_exclusive_unpaired_bases
                  << " Unband Count: " << node.number_of_children_outside_band
                  << " Children Count: " << node.children.size() << std::endl;
        energy += Gw;
        energy += Gwh * (node.number_of_bands - 2);
        energy += P_tilda * 2 * node.number_of_bands;
        energy += Q_tilda * 2 * node.number_of_exclusive_unpaired_bases;
        energy += P_i * node.number_of_children_outside_band;
        energy *= 100;

        return energy;
    }

   private:
};
}  // namespace knotergy