#pragma once

#include "../loops/LoopNode.hpp"

namespace knotergy {
class PseudoknotFunctions {
   public:
    PseudoknotFunctions();
    ~PseudoknotFunctions();

   private:
    int pseudoknot_energy(LoopNode node) {
        float P_tilda = 0.1;  // pair in a pseudoknot
        float Q_tilda = 0.2;  // unpaired bases in a pseudoknot
        float P_i = 0.1;      // for E&R energy model
        float Gw = 7;         // starting a pseudoknot loop
        float Gwh = 6;        // each more band region

        float energy;
    }
};
}  // namespace knotergy