#include "BandFinder.hpp"

namespace knotergy {

// Returns all bands within the specified region
std::vector<Band> BandFinder::find_bands(size_t cr_start, size_t cr_end, LoopType loop_type,
                                         std::vector<PairedBaseNode>& aux_bands,
                                         const std::vector<size_t>& pair_table,
                                         const std::vector<size_t>& cr_pair_table) {
    // sanity check bounds
    if (cr_end >= pair_table.size()) THROW_ERROR("Right bound exceeds the size of structure.");

    // If not a pseudoknot, there are no bands to find (bands only exist in pseudoknots)
    if (loop_type != LoopType::Pseudoknot) {
        return {};
    }

    std::vector<Band> bands;
    bands.reserve(8);  // Most pseudoknots are small, so reserving some just in case

    // Linked list pointing to potential band boundaries
    generate_paired_base_links(cr_start, cr_end, aux_bands, pair_table, cr_pair_table);

    // Iterate through linked paired bases only
    size_t band_start = cr_start;

    while (band_start < cr_end) {
        // Skips closing base pairs (Bands start on the opening base of a pair)
        if (pair_table[band_start] < band_start) {
            band_start = aux_bands[band_start].next;
            continue;
        }

        // pairing outside of closed region (means invalid input, region isn't really closed)
        if (pair_table[band_start] > cr_end) {
            THROW_ERROR("Pairing outside of \"closed region\" detected." +
                        std::to_string(band_start) + " pairs with " +
                        std::to_string(pair_table[band_start]) + " which is outside [" +
                        std::to_string(cr_start) + ", " + std::to_string(cr_end) + "].");
        }

        size_t band_end = pair_table[band_start];

        // Finds the inner band positions
        auto [left_inner, right_inner] =
            find_stem_inner_indices(band_start, band_end, aux_bands, pair_table);

        // Create the band and add to list
        bands.emplace_back(band_start, left_inner, right_inner, band_end, pair_table,
                           cr_pair_table);

        band_start = aux_bands[left_inner].next;
    }

    return bands;
}
// Convenience method for LoopNode
std::vector<Band> BandFinder::find_bands(const LoopNode& node,
                                         std::vector<PairedBaseNode>& aux_bands,
                                         const ProcessedRNAEntry& processed_rna) {
    return find_bands(node.begin, node.end, node.loop_type, aux_bands,
                      processed_rna.get_pair_table(),
                      processed_rna.get_closed_regions_pair_table());
}

// Extends the stem to find inner band positions
std::pair<size_t, size_t> BandFinder::find_stem_inner_indices(
    size_t band_start, size_t band_end, const std::vector<PairedBaseNode>& aux_bands,
    const std::vector<size_t>& pair_table) {
    size_t left_inner = band_start;
    size_t right_inner = band_end;

    while (true) {
        size_t next_left_inner = aux_bands[left_inner].next;
        size_t next_right_inner = aux_bands[right_inner].prev;

        if (next_left_inner == NULL_INDEX || next_right_inner == NULL_INDEX) {
            THROW_ERROR("PairedBaseNode points to NULL_INDEX.");
        }

        // If the next inner positions cross, we have reached the end of the stem
        if (next_left_inner >= next_right_inner) {
            break;
        }

        // If the next inner positions are not paired with each other, we have reached the end of
        // the stem
        if (pair_table[next_left_inner] != next_right_inner) {
            break;
        }

        left_inner = next_left_inner;
        right_inner = next_right_inner;
    }

    return {left_inner, right_inner};
}

void BandFinder::generate_paired_base_links(size_t cr_start, size_t cr_end,
                                            std::vector<PairedBaseNode>& aux_bands,
                                            const std::vector<size_t>& pair_table,
                                            const std::vector<size_t>& cr_pair_table) {
    if (cr_end < cr_start) return;

    if (pair_table[cr_start] == NULL_INDEX) {
        THROW_ERROR("Left index is not a base-pair");
    }
    if (pair_table[cr_end] == NULL_INDEX) {
        THROW_ERROR("Right index is not a base-pair");
    }
    if (cr_pair_table[cr_start] != cr_end) {
        THROW_ERROR("Left bound and right bound do not form a closed region");
    }

    size_t prev_key = cr_start;
    for (size_t i = cr_start + 1; i < cr_end; ++i) {
        // skip unpaired
        if (pair_table[i] == NULL_INDEX) {
            continue;
        }

        // skip nested closed regions
        if (cr_pair_table[i] != NULL_INDEX) {
            i = cr_pair_table[i];
            continue;
        }

        aux_bands[prev_key].next = i;
        aux_bands[i].prev = prev_key;
        prev_key = i;
    }

    // link the last one to the end of the closed region
    aux_bands[prev_key].next = cr_end;
    aux_bands[cr_end].prev = prev_key;
    return;
}

void BandFinder::generate_paired_base_links(const LoopNode& node,
                                            std::vector<PairedBaseNode>& aux_bands,
                                            const ProcessedRNAEntry& processed_entry) {
    generate_paired_base_links(node.begin, node.end, aux_bands, processed_entry.get_pair_table(),
                               processed_entry.get_closed_regions_pair_table());
}

}  // namespace knotergy