#pragma once

#include <vector>

#include "RNAEntry.hpp"

namespace ComputeEnergy {

/**
 * @brief Metadata for a single position in the RNA band.
 *
 * Represents a node in a doubly linked list used to track band regions
 * within a pseudoknotted RNA structure.
 *
 * - `is_band_start` marks the left (starting) boundary of a band.
 * - `band_end` points to the right boundary of the band.
 * - `next` and `prev` allow traversal through the candidate positions.
 */

struct B_pattern {
    bool is_band_start;  // checks if it's the start of a band region (opening bracket)
    size_t next, prev;   // next and prev index in pattern in RNAEntry
    size_t band_end;     // right border of the band region (closing bracket)
};

class Bands {
   public:
    Bands(RNAEntry& entry);

    void update_links(size_t from, size_t to);
    size_t find_next_good_index(size_t start_index, size_t end_index);

    void aux_Find_bands(size_t border1, size_t border2, size_t* NumberOfBands,
                        size_t* CurrentBandRegion);
    size_t prev(size_t i);
    size_t next(size_t i);
    std::vector<B_pattern> pattern;
    RNAEntry entry;
};

}  // namespace ComputeEnergy