#include "BandFinder.hpp"

namespace knotergy {

// Returns all bands within the specified region
// Returns all bands within the specified region
std::vector<Band> BandFinder::find_bands(size_t left_bound, size_t right_bound,
                                         LoopType loop_type, std::vector<PairedBaseNode>& aux_bands,
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
    generate_paired_base_links(left_bound, right_bound, aux_bands, pairings, cr_pairings);

    // Iterate through linked paired bases only
    size_t i = left_bound;

    while (i <= right_bound) {
        // Skips closing base pairs
        if (pairings[i] < i) {
            i = aux_bands[i].next;
            continue;
        }

        // pairing outside of closed region (means invalid input, region isn't really closed)
        if (pairings[i] > right_bound) {
            THROW_ERROR("Pairing outside of \"closed region\" detected." + std::to_string(i) +
                        " pairs with " + std::to_string(pairings[i]) + " which is outside [" +
                        std::to_string(left_bound) + ", " + std::to_string(right_bound) + "].");
        }

        size_t j = pairings[i]; 
        size_t i_prime = i;      // will be left inner after extension
        size_t j_prime = j;      // will be right inner after extension

        // Traverses the band until the innermost band positions (i`, j`) are found
        while (extend_stem(i_prime, j_prime, aux_bands, pairings)) {}

        bands.push_back(Band{i, i_prime, j_prime, j, pairings, cr_pairings});

        i = aux_bands[i_prime].next;
    }

    return bands;
}
// Convenience method for LoopNode
std::vector<Band> BandFinder::find_bands(const LoopNode& node,
                                         std::vector<PairedBaseNode>& aux_bands,
                                         const ProcessedRNAEntry& processed_rna) {
    return find_bands(node.begin, node.end, node.loop_type, aux_bands, processed_rna.get_pairings(),
                      processed_rna.get_closed_regions_pairings());
}

// Extends the stem to find inner band positions
bool BandFinder::extend_stem(size_t& i_prime, size_t& j_prime, std::vector<PairedBaseNode>& aux_bands,
                             const std::vector<size_t>& pairings) {
    
    // Jump to next paired base
    size_t i_tmp = aux_bands[i_prime].next;
    size_t j_tmp = aux_bands[j_prime].prev;

    // i should be before j (terminate before they cross)
    if (i_tmp >= j_tmp) {
        return false;
    }

    // sanity check (should never happen)
    if (i_tmp == NULL_INDEX || j_tmp == NULL_INDEX) {
        THROW_ERROR("PairedBaseNode points to NULL_INDEX.");
    }

    if (pairings[i_tmp] != j_tmp) {
        return false;
    }

    i_prime = i_tmp;
    j_prime = j_tmp;
    return true;
}

void BandFinder::generate_paired_base_links(
    size_t left_bound,
    size_t right_bound,
    std::vector<PairedBaseNode>& aux_bands,
    const std::vector<size_t>& pairings,
    const std::vector<size_t>& cr_pairings) 
{
    if (right_bound < left_bound) return;

    if (pairings[left_bound] == NULL_INDEX) {
        THROW_ERROR("Left index is not a base-pair");
    }
    if (pairings[right_bound] == NULL_INDEX) {
        THROW_ERROR("Right index is not a base-pair");
    }
    if (cr_pairings[left_bound] != right_bound) {
        THROW_ERROR("Left bound and right bound do not form a closed region");
    }
    
    size_t prev_key = left_bound;
    for (size_t i = left_bound + 1; i < right_bound; ++i) {

        // skip unpaired
        if (pairings[i] == NULL_INDEX) {
            continue;
        }

        // skip nested closed regions
        if (cr_pairings[i] != NULL_INDEX) {
            i = cr_pairings[i];
            continue;
        }

        aux_bands[prev_key].next = i;
        aux_bands[i].prev = prev_key;
        prev_key = i;
    }

    aux_bands[prev_key].next = right_bound;
    aux_bands[right_bound].prev = prev_key;
    return;
}

void BandFinder::generate_paired_base_links(
    const LoopNode& node, std::vector<PairedBaseNode>& aux_bands,
    const ProcessedRNAEntry& processed_entry) {
    generate_paired_base_links(node.begin, node.end, aux_bands, processed_entry.get_pairings(),
                                      processed_entry.get_closed_regions_pairings());
}

}  // namespace knotergy