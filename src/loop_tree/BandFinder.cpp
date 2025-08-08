#include "BandFinder.hpp"

namespace knotergy {



std::vector<Band> BandFinder::find_bands(const size_t& left_bound, const size_t& right_bound,
                                         const LoopType& loop_type,
                                         const std::vector<size_t>& pairings,
                                         const std::vector<size_t>& cr_pairings) {
    std::vector<Band> bands;

    // std::unordered_map<size_t, size_t> next_valid_bp = generate_next_bp_map(left_bound, right_bound, pairings, cr_pairings);
    // std::unordered_map<size_t, size_t> prev_valid_bp = generate_prev_bp_map(left_bound, right_bound, pairings, cr_pairings);

    if (loop_type != LoopType::Pseudoknot) {
        bands.push_back(Band{left_bound, left_bound, right_bound, right_bound, pairings});
        return bands;
    }

    const size_t n = pairings.size();
    if (right_bound >= n) {
        THROW_ERROR("Right bound exceeds the size of structure.");
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
        while (extend_stem(i_prime, j_prime, left_bound, right_bound, pairings)) {
        }

        // stores entire band
        bands.push_back(Band{i, i_prime, j_prime, j, pairings});

        i = i_prime;  // fast-forward
    }
    return bands;
}

std::vector<Band> BandFinder::find_bands(const LoopNode& node, const RNAProcessedEntry& processed_rna) {
    return find_bands(node.begin, node.end, node.loop_type, processed_rna.get_pairings(), processed_rna.get_closed_regions_pairings());
}

bool BandFinder::extend_stem(size_t& i_prime, size_t& j_prime, const size_t& left_bound,
                             const size_t& right_bound, const std::vector<size_t>& pairings) {
    size_t i_tmp = i_prime + 1;
    std::cout << "i: " << i_tmp << std::endl;
    size_t j_tmp = j_prime - 1;
    std::cout << "j: " << j_tmp << std::endl;

    /* skip unpaired positions on either side */
    while (i_tmp <= right_bound && pairings[i_tmp] == NULL_INDEX) {
        std::cout << "Skip Forward: " << i_tmp << std::endl;
        ++i_tmp;
    }
    while (j_tmp > left_bound && pairings[j_tmp] == NULL_INDEX) {
        std::cout << "Skip Backwards: " << i_tmp << std::endl;
        --j_tmp;
    }

    if (i_tmp < j_tmp && pairings[i_tmp] == j_tmp) {  // still a canonical stack
        i_prime = i_tmp;
        j_prime = j_tmp;
        std::cout << "True" << std::endl;
        return true;
    }
    std::cout << "False" << std::endl;
    return false;
}

std::pair<std::shared_ptr<LinkedList>, std::shared_ptr<LinkedList>> const BandFinder::generate_aux_band_list(const size_t& left_bound,
                                                                          const size_t& right_bound,
                                                                          const std::vector<size_t>& pairings,
                                                                          const std::vector<size_t>& cr_pairings){
    if (right_bound < left_bound) return std::make_pair(nullptr, nullptr);

    size_t start_idx = left_bound;

    while (pairings[start_idx] == NULL_INDEX){
        ++start_idx;
    }
    
    std::shared_ptr<LinkedList> start = std::make_shared<LinkedList>(start_idx);
    std::shared_ptr<LinkedList> current = start;

    for (size_t i = start_idx + 1; i <= right_bound; ++i){

        if (pairings[i] != NULL_INDEX){
            continue;
        }

        if (cr_pairings[i] != NULL_INDEX && i != left_bound) {
            i = cr_pairings[i];
            continue;
        }

        std::shared_ptr<LinkedList> new_base = std::make_shared<LinkedList>(i);
        new_base->prev = current;
        current->next = new_base;
        current = new_base;
    }

    std::shared_ptr<LinkedList> end = current;
    return std::make_pair(start, end);
}
std::pair<std::shared_ptr<LinkedList>, std::shared_ptr<LinkedList>> const BandFinder::generate_aux_band_list(const LoopNode& node, const RNAProcessedEntry& processed_entry){
    return generate_aux_band_list(node.begin, node.end, processed_entry.get_pairings(), processed_entry.get_closed_regions_pairings());
}

}  // namespace knotergy