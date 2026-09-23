#pragma once

#include <ViennaRNA/params/basic.hpp>

namespace knotergy {
[[nodiscard]] inline constexpr int add_or_inf(int a, int b) {
    if (a >= INF || b >= INF) {
        return INF;
    }
    return a + b;
}

/**
 * @brief Stores dangle energy values for all four dangle configurations.
 *
 * Represents the energy contributions of dangling ends (unpaired nucleotides adjacent
 * to base pairs) in various configurations.
 */
struct DangleSet {
   public:
    DangleSet() : no_dangle(0), left_dangle(0), right_dangle(0), both_dangle(0) {}
    DangleSet(int no_dangle_energy, int left_dangle_energy, int right_dangle_energy,
              int both_dangle_energy)
        : no_dangle(no_dangle_energy),
          left_dangle(left_dangle_energy),
          right_dangle(right_dangle_energy),
          both_dangle(both_dangle_energy) {}

    int no_dangle;     ///< Energy with no dangling ends.
    int left_dangle;   ///< Energy with only 5' dangling end.
    int right_dangle;  ///< Energy with only 3' dangling end.
    int both_dangle;   ///< Energy with both dangling ends.

    /**
     * @brief Get the minimum (most favorable) energy among all dangle configurations.
     *
     * @return The minimum energy value in centicalories.
     */
    [[nodiscard]] int best() const {
        return std::min({no_dangle, left_dangle, right_dangle, both_dangle});
    }

    /**
     * @brief Get the minimum energy between no dangle and left dangle only.
     *
     * @return The minimum energy value in centicalories.
     */
    [[nodiscard]] int best_left() const { return std::min(no_dangle, left_dangle); }

    /**
     * @brief Get the minimum energy between no dangle and right dangle only.
     *
     * @return The minimum energy value in centicalories.
     */
    [[nodiscard]] int best_right() const { return std::min(no_dangle, right_dangle); }

    /** @brief Add a constant energy to all dangle configurations.
     *
     * This is useful for applying energy contributions that affect all configurations
     * equally, such as terminal mismatches
     *
     * @param energy The energy value in centicalories to add to each configuration.
     * @return Reference to the modified DangleSet.
     */
    DangleSet& operator+=(int energy) {
        no_dangle += energy;
        left_dangle += energy;
        right_dangle += energy;
        both_dangle += energy;
        return *this;
    }
};

}  // namespace knotergy