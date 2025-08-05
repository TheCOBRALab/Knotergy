#pragma once

#include <cstddef>
#include <vector>

#include "../loop_tree/LoopNode.hpp"
#include "../pipeline/shared.hpp"
#include "Band.hpp"
namespace knotergy {

// walk one step along a perfect stack:  (i+1) pairs (j−1)
class BandFinder {
   public:
    BandFinder(const std::vector<size_t>& pairings, const std::vector<size_t>& closed_region_pairings)
        : pairings_{pairings},
          closed_region_pairings_{closed_region_pairings},
          done_(pairings.size(), false) {};

    BandFinder(const RNAProcessedEntry& processed_rna)
        : BandFinder(processed_rna.get_pairings(), processed_rna.get_closed_region_pairings()) {};
    
    /*──────────────── attach bands + pseudo-nest info to one node ────────────*/
    void annotate_bands(LoopNode& node) {
        node.bands = find_bands_in_region(node);
        node.number_of_bands = static_cast<int>(node.bands.size());

        /* only pseudoknots need bands; leave the rest untouched */
        if (node.loop_type != LoopType::Pseudoknot) return;

        const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
        const std::vector<Band>& bands = node.bands;
        size_t child_idx = 0;
        size_t band_idx = 0;
        
        // checks if a child is exclusively in one band (InsideBand), or if its CrossBand 
        while (child_idx < children.size()) {
            std::shared_ptr<LoopNode> child = children[child_idx];
            child->pseudo_type = PseudoNestedType::CrossBand;

            // Prevent out-of-bounds in band lookup
            if (band_idx >= bands.size()) {
                ++child_idx;
                continue;
            }

            // Checks if child is exclusively in one band
            if (node.bands[band_idx].contains(child->begin, child->end)) {
                bool crosses_previous = 
                (band_idx > 0) && bands[band_idx - 1].contains(child->begin, child->end);
 
                bool crosses_next = 
                (band_idx + 1) < bands.size() && bands[band_idx + 1].contains(child->begin, child->end);

                if (!crosses_previous && !crosses_next) {
                    child->pseudo_type = PseudoNestedType::InsideBand;
                }
                ++child_idx;  // Done with this child, move to next
            } else {
                ++band_idx;
            }
        }
    }

    void annotate_bands(const std::shared_ptr<LoopNode>& node){
        annotate_bands(*node);
    }


   private:
    const std::vector<size_t>& pairings_;
    const std::vector<size_t>& closed_region_pairings_;
    std::vector<bool> done_;

bool extend_stem(size_t& i_prime, size_t& j_prime, const size_t& left_bound,
                        const size_t& right_bound) {
            size_t i_tmp = i_prime + 1;
            std::cout << "i: " << i_tmp << std::endl;
            size_t j_tmp = j_prime - 1;
            std::cout << "j: " << j_tmp << std::endl;

            /* skip unpaired positions on either side */
            while (i_tmp <= right_bound && pairings_[i_tmp] == NULL_INDEX) {
                std::cout << "Skip Forward: " << i_tmp << std::endl;
                ++i_tmp;
            }
            while (j_tmp > left_bound && pairings_[j_tmp] == NULL_INDEX) {
                std::cout << "Skip Backwards: " << i_tmp << std::endl;
                --j_tmp;
            }

            if (i_tmp < j_tmp && pairings_[i_tmp] == j_tmp) {  // still a canonical stack
                i_prime = i_tmp;
                j_prime = j_tmp;
                std::cout << "True" << std::endl;
                return true;
            }
            std::cout << "False" << std::endl;
            return false;
        }

    /*──────────────── find all bands inside [left, right] ────────────────────*/
    std::vector<Band> find_bands_in_region(LoopNode node) {
        size_t left_bound = node.begin;
        size_t right_bound = node.end;
        std::vector<Band> bands;
        

        if (node.loop_type != LoopType::Pseudoknot){
            done_[right_bound] = true;
            done_[left_bound] = true;
            bands.push_back(Band{left_bound, left_bound, right_bound, right_bound, pairings_});
            return bands;
        }

        const size_t n = pairings_.size();
        if (right_bound >= n) {
            right_bound = n - 1;
        }

        for (size_t i = left_bound; i <= right_bound; ++i) {
            /* skip anything that is   – already used
             *                         – unpaired
             *                         – closing half
             *                         – pairs outside this region */
            if (pairings_[i] == NULL_INDEX || pairings_[i] < i ||
                pairings_[i] > right_bound)
                continue;
            
            if (closed_region_pairings_[i] != NULL_INDEX && i != left_bound){
                i = closed_region_pairings_[i];
                continue;
            }

            size_t j = pairings_[i];
            size_t i_prime = i;
            size_t j_prime = j;

            // walks the stem until last base pair in the given region (finds i` and j`)
            while (extend_stem(i_prime, j_prime, left_bound, right_bound)) {}

            // stores entire band
            bands.push_back(Band{i, i_prime, j_prime, j, pairings_});

            /* mark every position that belongs to this band */
            for (size_t k = i; k <= i_prime; ++k) {
                done_[k] = true;
                if (pairings_[k] != NULL_INDEX) done_[pairings_[k]] = true;
                if (closed_region_pairings_[k]!= NULL_INDEX && closed_region_pairings_[k] > k){
                    k = closed_region_pairings_[k];
                }
            }
            for (size_t k = j_prime; k <= j; ++k) {
                done_[k] = true;
                if (pairings_[k] != NULL_INDEX) done_[pairings_[k]] = true;
                if (closed_region_pairings_[k]!= NULL_INDEX && closed_region_pairings_[k] > k){
                    k = closed_region_pairings_[k];
                }
            }

            i = i_prime;  // fast-forward
        }
        return bands;
    }

    // std::vector<Band> find_bands_in_region(size_t start) {
    //     return find_bands_in_region(start, pairings_.size() - 1);
    // }

    // std::vector<Band> find_bands_in_region() {
    //     return find_bands_in_region(0, pairings_.size() - 1);
    // }
};

}  // namespace knotergy