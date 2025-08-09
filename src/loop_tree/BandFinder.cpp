#include "BandFinder.hpp"

namespace knotergy {



std::vector<Band> BandFinder::find_bands(const size_t& left_bound, const size_t& right_bound,
                                         const LoopType& loop_type,
                                         const std::vector<size_t>& pairings,
                                         const std::vector<size_t>& cr_pairings) {
    std::vector<Band> bands;

    if (loop_type != LoopType::Pseudoknot) {
        bands.push_back(Band{left_bound, left_bound, right_bound, right_bound, pairings});
        return bands;
    }

    const size_t n = pairings.size();
    if (right_bound >= n) {
        THROW_ERROR("Right bound exceeds the size of structure.");
    }

    std::unordered_map<size_t, Aux_Band> aux_bands = generate_Aux_Band_list(left_bound, right_bound, pairings, cr_pairings);
    
    if (aux_bands.empty()) {
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
        while (extend_stem(i_prime, j_prime, aux_bands, pairings)) {
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

bool BandFinder::extend_stem(size_t& i_prime, size_t& j_prime, const std::unordered_map<size_t, Aux_Band>& aux_bands, const std::vector<size_t>& pairings) {
    
    if (aux_bands.find(i_prime) != aux_bands.end() && aux_bands.find(j_prime) != aux_bands.end()) {
        size_t i_tmp = aux_bands.at(i_prime).next;
        size_t j_tmp = aux_bands.at(j_prime).prev;
        std::cout << "i_tmp: " << i_tmp << std::endl;
        std::cout << "j_tmp: " << j_tmp << std::endl;

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

std::unordered_map<size_t, Aux_Band> const BandFinder::generate_Aux_Band_list(const size_t& left_bound,
                                                                          const size_t& right_bound,
                                                                          const std::vector<size_t>& pairings,
                                                                          const std::vector<size_t>& cr_pairings){
    if (right_bound < left_bound) return {};
    

    
    size_t start_idx = left_bound;

    // finds the first paired base in the region (Normally it's paired, but just in case)
    [[unlikely]] while (pairings[start_idx] == NULL_INDEX){
        ++start_idx;
    }
    
    std::unordered_map<size_t, Aux_Band> aux_bands;
    aux_bands.emplace(start_idx, Aux_Band(start_idx));
    size_t current_key = start_idx;

    for (size_t i = start_idx + 1; i < right_bound; ++i){

        if (pairings[i] == NULL_INDEX){
            continue;
        }

        if (cr_pairings[i] != NULL_INDEX && i != left_bound) {
            i = cr_pairings[i];
            continue;
        }

        if (aux_bands.find(current_key) == aux_bands.end()) {
            THROW_ERROR("Current key not found in aux_bands");
        }
        
        std::cout << "Adding aux band for i: " << i << std::endl;
        aux_bands[current_key].next = i;
        aux_bands[i] = Aux_Band(i, current_key);
        current_key = i;
    }

    // Normally the right_bound is paired, but just in case
    [[likely]] if (pairings[right_bound] != NULL_INDEX) {
        aux_bands[current_key].next = right_bound;
        aux_bands[right_bound] = Aux_Band(right_bound, current_key);
    }

    return aux_bands;
}
std::unordered_map<size_t, Aux_Band> const BandFinder::generate_Aux_Band_list(const LoopNode& node, const RNAProcessedEntry& processed_entry){
    return generate_Aux_Band_list(node.begin, node.end, processed_entry.get_pairings(), processed_entry.get_closed_regions_pairings());
}

}  // namespace knotergy