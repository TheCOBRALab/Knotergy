#include "Bands.hpp"

#include <vector>

#include "RNAEntry.hpp"

namespace ComputeEnergy {

Bands::Bands(RNAEntry& entry) : entry(entry) {
    for (size_t i = 1; i < entry.structure.size(); ++i) {
        pattern[i].is_band_start = false;
        pattern[i].prev = i - 1;
        pattern[i].next = i + 1;
        pattern[i].band_end = 0;
    }
    for (size_t i = 1; i < entry.structure.size(); ++i)
        if (entry.get_pairings()[i] <= 0) {
            update_links(pattern[i].prev, pattern[i].next);
        }
};

/*********************************************************************************
Find_next_good_index: returns the start of the next band region
*********************************************************************************/
size_t Bands::find_next_good_index(size_t start_index, size_t end_index) {
    size_t i = start_index;
    while ((pattern[i].is_band_start == true) && (i < end_index + 1)) {
        i = pattern[i].next;
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

void Bands::aux_Find_bands(size_t border1, size_t border2, size_t* NumberOfBands,
                           size_t* CurrentBandRegion) {
    size_t i, j, i_prime, j_prime;  // the band will be [i,i_prime] [j_prime,j]
    size_t i_help, j_help;
    i = find_next_good_index(border1, border2);
    *NumberOfBands = 0;
    *CurrentBandRegion = 0;
    while (i < border2 + 1) {
        *NumberOfBands = *NumberOfBands + 1;
        j = entry.get_pairings()[i];
        i_prime = i;
        j_prime = j;
        i_help = pattern[i_prime].next;
        j_help = pattern[j_prime].prev;
        while (entry.get_pairings()[i_help] == j_help) {
            i_prime = i_help;
            j_prime = j_help;
            i_help = pattern[i_prime].next;
            j_help = pattern[j_prime].prev;
        }

        // found the borders of the band region, now set them in pattern
        pattern[i].is_band_start = true;
        pattern[j_prime].is_band_start = true;
        pattern[i].band_end = i_prime;
        pattern[j_prime].band_end = j;
        if (*CurrentBandRegion < j_prime) *CurrentBandRegion = j_prime;
        if (i != i_prime) {  // if the band is not just a base pair
            update_links(
                i, pattern[i_prime]
                       .next);  // the next possible left border of a region of a band, after i, is
                                // now after i_prime, since [i,i_prime] is now known to be a band
            update_links(j_prime, pattern[j].next);  // the next possible left border of a region of
                                                     // a band, after j_prime, is now after j, since
                                                     // [j_prime,j] is now known to be a band
        }
        // check for the next possible basepairs in a band starting from i.next = i_prime.next (i.e.
        // starting from the next possible base pair outside of [i,i_prime])
        i = find_next_good_index(pattern[i].next, border2);
        fflush(stdout);
    }
}

void Bands::update_links(size_t from, size_t to) {
    pattern[from].next = to;
    pattern[to].prev = from;
};

size_t Bands::prev(size_t i) {
    return pattern[i].prev;
};

size_t Bands::next(size_t i) {
    return pattern[i].next;
};
}  // namespace ComputeEnergy