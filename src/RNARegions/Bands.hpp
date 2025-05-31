#pragma once

#include <iostream>
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
    bool is_band_start;  // if this index is the *left* border of a band region (either i or j')
    size_t next, prev;   // next and prev index in `pattern` in RNAEntry of the next pair
    size_t band_end;     // index of *right* end of this band region (i.e., i' for i, or j for j')

    void print_band() {
        std::cout << "is_band_start: " << is_band_start << std::endl
                  << "next: " << next << std::endl
                  << "prev: " << prev << std::endl
                  << "band_end: " << band_end << std::endl;
    }
};

class Bands {
   public:
    Bands(RNAEntry& entry_);

    void update_links(size_t from, size_t to);
    size_t find_next_good_index(size_t start_index, size_t end_index);

    std::pair<size_t, size_t> aux_find_bands(size_t border1, size_t border2);
    size_t prev(size_t i) const;
    size_t next(size_t i) const;
    void print_pairings();
    std::vector<B_pattern> pattern;
    RNAEntry entry;
};

}  // namespace ComputeEnergy