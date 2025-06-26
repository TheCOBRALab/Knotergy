#pragma once

namespace compute_energy {

class ViennaFunctions {
   public:
    static float stack_energy(int i, int j, int* sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    static float hairpin_energy(int i, int j, int* sequence, char* csequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    static float internal_loop_energy(int i, int j, int ip, int jp, int* sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    static float multi_energy(int i, int* sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    static float pseudoknot_energy(int i, int j, int* sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    static float external_energy(int i, int j, int* sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }
};
}  // namespace compute_energy