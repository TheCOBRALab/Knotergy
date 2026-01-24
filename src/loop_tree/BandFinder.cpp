#include "BandFinder.hpp"

namespace knotergy {

std::vector<Band> BandFinder::find_bands(const size_t& left_bound, const size_t& right_bound,
                                         const LoopType& loop_type,
                                         const std::vector<size_t>& pairings,
                                         const std::vector<size_t>& cr_pairings) {
    std::vector<Band> bands;

    if (loop_type != LoopType::Pseudoknot) {
        bands.push_back(
            Band{left_bound, left_bound, right_bound, right_bound, pairings, cr_pairings});
        return bands;
    }

    const size_t n = pairings.size();
    if (right_bound >= n) {
        THROW_ERROR("Right bound exceeds the size of structure.");
    }

    std::unordered_map<size_t, BandLink> band_links =
        generate_band_links(left_bound, right_bound, pairings, cr_pairings);

    if (band_links.empty()) {
        return bands;  // No bands found
    }

    for (size_t i = left_bound; i <= right_bound; ++i) {
        /* skip anything that is   – unpaired
         *                         – closing half
         *                         – pairs outside this region */
        if (pairings[i] == NULL_INDEX || pairings[i] < i || pairings[i] > right_bound) continue;

        // if nested closed region, skip
        if (cr_pairings[i] != NULL_INDEX && i != left_bound) {
            i = cr_pairings[i];
            continue;
        }

        size_t j = pairings[i];
        size_t i_prime = i;
        size_t j_prime = j;

        // walks the stem until last base pair in the given region (finds i` and j`)
        while (extend_stem(i_prime, j_prime, band_links, pairings)) {
        }

        // stores entire band
        bands.push_back(Band{i, i_prime, j_prime, j, pairings, cr_pairings});

        i = i_prime;  // fast-forward
    }
    return bands;
}

std::vector<Band> BandFinder::find_bands(const LoopNode& node,
                                         const ProcessedRNAEntry& processed_rna) {
    return find_bands(node.begin, node.end, node.loop_type, processed_rna.get_pairings(),
                      processed_rna.get_closed_regions_pairings());
}

bool BandFinder::extend_stem(size_t& i_prime, size_t& j_prime,
                             const std::unordered_map<size_t, BandLink>& band_links,
                             const std::vector<size_t>& pairings) {
    if (band_links.find(i_prime) != band_links.end() &&
        band_links.find(j_prime) != band_links.end()) {
        size_t i_tmp = band_links.at(i_prime).next;
        size_t j_tmp = band_links.at(j_prime).prev;

        if (i_tmp == NULL_INDEX || i_tmp >= j_tmp) {
            return false;  // no more extension possible
        }

        if (pairings[i_tmp] != j_tmp) {
            return false;  // not a valid base pair
        }

        i_prime = i_tmp;
        j_prime = j_tmp;
        return true;  // extension successful
    } else {
        return false;  // no aux band found for i or j
    }
}

std::unordered_map<size_t, BandLink> const BandFinder::generate_band_links(
    const size_t& left_bound, const size_t& right_bound, const std::vector<size_t>& pairings,
    const std::vector<size_t>& cr_pairings) {
    if (right_bound < left_bound) return {};

    if (pairings[left_bound] == NULL_INDEX) {
        THROW_ERROR("Left index is not a base-pair");
    }
    if (pairings[right_bound] == NULL_INDEX) {
        THROW_ERROR("Right index is not a base-pair");
    }
    if (cr_pairings[left_bound] != right_bound) {
        THROW_ERROR("Left bound does not pair with right bound");
    }

    std::unordered_map<size_t, BandLink> band_links;
    band_links.emplace(left_bound, BandLink{left_bound});
    size_t current_key = left_bound;

    for (size_t i = left_bound + 1; i < right_bound; ++i) {
        // skip unpaired
        if (pairings[i] == NULL_INDEX) {
            continue;
        }

        // skip closed regions
        if (cr_pairings[i] != NULL_INDEX && i != left_bound) {
            i = cr_pairings[i];
            continue;
        }

        if (band_links.find(current_key) == band_links.end()) {
            THROW_ERROR("Current key not found in band_links");
        }

        band_links[current_key].next = i;
        band_links[i] = BandLink{i, current_key};
        current_key = i;
    }
    // process last link
    band_links[current_key].next = right_bound;
    band_links[right_bound] = BandLink{right_bound, current_key};

    return band_links;
}
std::unordered_map<size_t, BandLink> const BandFinder::generate_band_links(
    const LoopNode& node, const ProcessedRNAEntry& processed_entry) {
    return generate_band_links(node.begin, node.end, processed_entry.get_pairings(),
                               processed_entry.get_closed_regions_pairings());
}

}  // namespace knotergy