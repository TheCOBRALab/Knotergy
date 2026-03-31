#include "BandFinder.hpp"

namespace knotergy {

// Returns all bands within the specified region
std::vector<Band> BandFinder::find_bands(const size_t& left_bound, const size_t& right_bound,
                                         const LoopType& loop_type,
                                         const std::vector<size_t>& pairings,
                                         const std::vector<size_t>& cr_pairings) {
    // sanity check bounds
    const size_t n = pairings.size();
    if (right_bound >= n) THROW_ERROR("Right bound exceeds the size of structure.");

    // If not a pseudoknot, there are no bands to find (bands only exist in pseudoknots)
    if (loop_type != LoopType::Pseudoknot) {
        return {};
    }

    std::vector<Band> bands;

    // Linked list pointing to potential band boundaries
    std::unordered_map<size_t, PairedBaseNode> paired_base_links =
        generate_paired_base_links(left_bound, right_bound, pairings, cr_pairings);

    if (paired_base_links.empty()) {
        return bands;  // No paired bases found? (should never happen)
    }

    // Iterate through the region to find bands
    for (size_t i = left_bound; i <= right_bound; ++i) {
        /* skip anything that is   – unpaired
         *                         – closing base */
        if (pairings[i] == NULL_INDEX || pairings[i] < i) continue;

        // pairing outside of closed region (means invalid input, region isn't really closed)
        if (pairings[i] > right_bound) {
            THROW_ERROR("Pairing outside of \"closed region\" detected." + std::to_string(i) +
                        " pairs with " + std::to_string(pairings[i]) + " which is outside [" +
                        std::to_string(left_bound) + ", " + std::to_string(right_bound) + "].");
        }

        // if nested closed region, skip (They are not part of the pseudoknot)
        if (cr_pairings[i] != NULL_INDEX && i != left_bound) {
            i = cr_pairings[i];
            continue;
        }

        size_t j = pairings[i];  // j is the paired base of i
        size_t i_prime = i;      // will be left inner after extension
        size_t j_prime = j;      // will be right inner after extension

        // walks the stem until last base pair in the given region (finds i` and j`)
        // Finds the full band that i and j belong to
        while (extend_stem(i_prime, j_prime, paired_base_links, pairings)) {
        }

        // stores entire bands
        bands.push_back(Band{i, i_prime, j_prime, j, pairings, cr_pairings});

        i = i_prime;  // fast-forward to inner position (since everything in between is part of this
                      // band)
    }
    return bands;
}

// Convenience method for LoopNode
std::vector<Band> BandFinder::find_bands(const LoopNode& node,
                                         const ProcessedRNAEntry& processed_rna) {
    return find_bands(node.begin, node.end, node.loop_type, processed_rna.get_pairings(),
                      processed_rna.get_closed_regions_pairings());
}

// Extends the stem to find inner band positions
bool BandFinder::extend_stem(size_t& i_prime, size_t& j_prime,
                             const std::unordered_map<size_t, PairedBaseNode>& paired_base_links,
                             const std::vector<size_t>& pairings) {
    auto it_i = paired_base_links.find(i_prime);
    if (it_i == paired_base_links.end()) return false;  // (should never happen if input is valid)

    auto it_j = paired_base_links.find(j_prime);
    if (it_j == paired_base_links.end()) return false;  // (should never happen if input is valid)

    // Jump to next paired base
    size_t i_tmp = it_i->second.next;
    size_t j_tmp = it_j->second.prev;

    // i should be before j (terminate before they cross)
    if (i_tmp >= j_tmp) {
        return false;
    }

    // sanity check
    if (i_tmp == NULL_INDEX) {
        THROW_ERROR("PairedBaseNode points to NULL_INDEX. Index: " + std::to_string(i_prime));
    }

    // validate that they are indeed paired
    if (pairings[i_tmp] != j_tmp) {
        return false;  // not a valid base pair
    }

    i_prime = i_tmp;
    j_prime = j_tmp;
    return true;  // extension successful, continue extending
}

// Linked list of base pair positions (Skips unpaired and closed regions)
std::unordered_map<size_t, PairedBaseNode> const BandFinder::generate_paired_base_links(
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
        THROW_ERROR("Left bound and right bound do not form a closed region");
    }

    // create linked list of band positions
    std::unordered_map<size_t, PairedBaseNode> paired_base_links;
    paired_base_links.emplace(left_bound, PairedBaseNode{left_bound});
    size_t prev_key = left_bound;

    // Traverse the region to build links
    for (size_t i = left_bound + 1; i < right_bound; ++i) {
        // skip unpaired
        if (pairings[i] == NULL_INDEX) {
            continue;
        }

        // skip closed regions
        if (cr_pairings[i] != NULL_INDEX) {
            i = cr_pairings[i];
            continue;
        }

        auto it = paired_base_links.find(prev_key);
        if (it == paired_base_links.end()) {
            THROW_ERROR("Prev key not found in paired_base_links");
        }

        // link to next base-pair position
        it->second.next = i;
        paired_base_links[i] = PairedBaseNode{i, prev_key};
        prev_key = i;
    }

    // process last link
    paired_base_links[prev_key].next = right_bound;
    paired_base_links[right_bound] = PairedBaseNode{right_bound, prev_key};

    return paired_base_links;
}
std::unordered_map<size_t, PairedBaseNode> const BandFinder::generate_paired_base_links(
    const LoopNode& node, const ProcessedRNAEntry& processed_entry) {
    return generate_paired_base_links(node.begin, node.end, processed_entry.get_pairings(),
                                      processed_entry.get_closed_regions_pairings());
}

}  // namespace knotergy