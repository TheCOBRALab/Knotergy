#include "Bands.hpp"

#include <vector>

#include "RNAEntry.hpp"

namespace ComputeEnergy {

Bands::Bands(RNAEntry& entry_) : entry(entry_) {
    // Index 0 is unused to simplify boundary logic — acts like a dummy head node
    size_t n = entry.get_structure().size();
    pattern.resize(n + 1);
    pattern[0] = B_pattern{false, 0, 0, 0};

    for (size_t i = 1; i <= n; ++i) {
        pattern[i] = B_pattern{false, (i + 1 <= n ? i + 1 : i), i - 1, 0};
    }

    for (size_t i = 1; i <= n; ++i) {
        if (entry.get_pairings()[i - 1] < 0) {
            update_links(pattern[i].prev, pattern[i].next);
        }
    }
}

/*********************************************************************************
Find_next_good_index: returns the start of the next band region
*********************************************************************************/
size_t Bands::find_next_good_index(size_t start, size_t end) {
    const std::size_t max_idx = entry.get_structure().size();
    const std::size_t last_safe = std::min(end, max_idx - 1);

    size_t i = start;
    while (i <= last_safe && pattern[i].is_band_start) {
        i = pattern[i].next;  // Assumes `next` is always valid
    }
    return i;
}

/*********************************************************************************
aux_Find_bands: returns the band regions of the pseudoknotted closed region,
Number of Bands and the pointer to the last band region.
*********************************************************************************/

// A band is the union of two regions [i,i_prime] and [j_prime,j]
// i and j_prime are considered to be left borders of the band region (isLeftBorder true)
// i_prime is considered to be the RightBorder of pattern[i]
// j is considered to be the RightBorder of pattern[j_prime]

std::pair<size_t, size_t> Bands::aux_find_bands(size_t border1, size_t border2) {
    size_t i, j, i_prime, j_prime;  // the bands will be [i,i_prime] [j_prime,j]
    size_t next_left_candidate, next_right_candidate;
    size_t band_count{0}, current_band_region{0};

    const std::vector<int>& pairings = entry.get_pairings();
    i = find_next_good_index(border1, border2);

    while (i < border2 + 1) {
        // find the borders i, i_prime, j, j_prime
        ++band_count;
        if (pairings[i] < 0) {
            throw std::runtime_error("Error: unpaired base in \"good index\"");
        }
        j = static_cast<size_t>(pairings[i]);
        i_prime = i;
        j_prime = j;
        next_left_candidate = pattern[i_prime].next;
        next_right_candidate = pattern[j_prime].prev;
        while (pairings[next_left_candidate] == static_cast<int>(next_right_candidate)) {
            i_prime = next_left_candidate;
            j_prime = next_right_candidate;
            next_left_candidate = pattern[i_prime].next;
            next_right_candidate = pattern[j_prime].prev;
        }

        // found the borders of the band region, now set them in pattern
        pattern[i].is_band_start = true;
        pattern[j_prime].is_band_start = true;
        pattern[i].band_end = i_prime;
        pattern[j_prime].band_end = j;
        if (current_band_region < j_prime) current_band_region = j_prime;
        if (i != i_prime) {  // if the band is not just a base pair

            // the next possible left border of a region of a band, after i, is now after i_prime,
            // since [i,i_prime] is now known to be a band
            update_links(i, pattern[i_prime].next);

            // the next possible left border of a region of a band, after j_prime, is now after j,
            // since [j_prime,j] is now known to be a band
            update_links(j_prime, pattern[j].next);
        }
        // check for the next possible basepairs in a band starting from i.next = i_prime.next (i.e.
        // starting from the next possible base pair outside of [i,i_prime])
        i = find_next_good_index(pattern[i].next, border2);
    }
    return {band_count, current_band_region};
}

void Bands::update_links(size_t from, size_t to) {
    pattern[from].next = to;
    pattern[to].prev = from;
}

size_t Bands::prev(size_t i) const {
    return pattern[i].prev;
}

size_t Bands::next(size_t i) const {
    return pattern[i].next;
}

void Bands::print_pairings() {
    size_t last_index_seen = 0;
    size_t i = 1;
    while (i > last_index_seen) {
        std::cout << i << std::endl;
        pattern[i].print_band();
        std::cout << std::endl;
        last_index_seen = i;
        i = pattern[i].next;
    }
}
}  // namespace ComputeEnergy