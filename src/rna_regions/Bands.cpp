#include "Bands.hpp"

#include <vector>

#include "../helpers/common.hpp"
#include "RNAEntry.hpp"

namespace compute_energy {

Bands::Bands(RNAEntry& entry) : entry_(entry) {
    size_t n = entry_.get_structure().size();
    pattern.resize(n + 1);  // includes dummy head at index 0

    // Index 0 acts like a dummy head node and tells us the index of the first pairing
    pattern[0] = B_pattern{false, (n > 0 ? 1 : NULL_INDEX), NULL_INDEX, 0};

    for (size_t i = 1; i < n; ++i) {
        pattern[i] = B_pattern{false, i + 1, i - 1, 0};
    }
    pattern[n] = B_pattern{false, NULL_INDEX, n - 1, 0};  // last element has no "next"

    for (size_t i = 1; i <= n; ++i) {
        if (entry_.get_pairings()[i - 1] == NULL_INDEX) {
            unlink(i);
        }
    }
}

/*********************************************************************************
Find_next_good_index: returns the start of the next band region
*********************************************************************************/
size_t Bands::find_next_good_index(size_t start, size_t end) {
    const std::size_t max_idx = entry_.get_structure().size();
    const std::size_t last_safe = (max_idx > 0) ? std::min(end, max_idx - 1) : 0;

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

    const std::vector<size_t>& pairings = entry_.get_pairings();
    i = find_next_good_index(border1, border2);

    while (i < border2 + 1) {
        // find the borders i, i_prime, j, j_prime
        ++band_count;
        if (pairings[i] == NULL_INDEX) {
            throw std::runtime_error("Error: unpaired base in \"good index\"");
        }
        j = pairings[i];
        i_prime = i;
        j_prime = j;
        next_left_candidate = pattern[i_prime].next;
        next_right_candidate = pattern[j_prime].prev;
        while (pairings[next_left_candidate] == next_right_candidate) {
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

void Bands::unlink(size_t idx) {
    if (pattern[idx].prev != NULL_INDEX) pattern[pattern[idx].prev].next = pattern[idx].next;
    if (pattern[idx].next != NULL_INDEX) pattern[pattern[idx].next].prev = pattern[idx].prev;
}

size_t Bands::prev(size_t i) const {
    return pattern[i].prev;
}

size_t Bands::next(size_t i) const {
    return pattern[i].next;
}

void Bands::print_pairings() {
    size_t i = pattern[0].next;
    while (i != NULL_INDEX) {
        pattern[i].print_band(i);
        i = pattern[i].next;
    }
}

}  // namespace compute_energy