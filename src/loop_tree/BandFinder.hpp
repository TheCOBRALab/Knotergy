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
    BandFinder(const std::vector<size_t>& pairings)
        : pairings_{pairings}, done_(pairings.size(), false) {};

    /*──────────────── attach bands + pseudo-nest info to one node ────────────*/
    void annotate_bands(const std::shared_ptr<LoopNode>& node) {
        node->bands = find_bands_in_region(node->begin, node->end);
        node->number_of_bands = static_cast<int>(node->bands.size());

        /* only pseudoknots need bands; leave the rest untouched */
        if (node->loop_type != LoopType::Pseudoknot) return;

        size_t child_idx = 0;
        size_t band_idx = 0;
        const std::vector<Band>& bands = node->bands;

        while (child_idx < node->children.size()) {
            std::shared_ptr<LoopNode>& child = node->children[child_idx];
            child->pseudo_type = PseudoNestedType::OutsideBand;

            // Prevent out-of-bounds in band lookup
            if (band_idx >= bands.size()) {
                ++child_idx;
                continue;
            }

            // Checks if it's exclusively in one band
            if (node->bands[band_idx].nests(child->begin, child->end)) {
                bool crosses_previous =
                    (band_idx > 0) && bands[band_idx - 1].nests(child->begin, child->end);
                bool crosses_next = (band_idx + 1 < bands.size()) &&
                                    bands[band_idx + 1].nests(child->begin, child->end);

                if (!crosses_previous && !crosses_next) {
                    child->pseudo_type = PseudoNestedType::InsideBand;
                }
                ++child_idx;  // Done with this child, move to next
            } else {
                ++band_idx;
            }
        }
    }

   private:
    const std::vector<size_t> pairings_;
    std::vector<bool> done_;

    bool extend_stem(size_t& i_prime, size_t& j_prime, const size_t& left_bound,
                     const size_t& right_bound) {
        size_t i_tmp = i_prime + 1;
        size_t j_tmp = j_prime - 1;

        /* skip unpaired positions on either side */
        while (i_tmp <= right_bound && pairings_[i_tmp] == NULL_INDEX) ++i_tmp;
        while (j_tmp > left_bound && pairings_[j_tmp] == NULL_INDEX) --j_tmp;

        if (i_tmp < j_tmp && pairings_[i_tmp] == j_tmp) {  // still a canonical stack
            i_prime = i_tmp;
            j_prime = j_tmp;
            return true;
        }
        return false;
    }

    /*──────────────── find all bands inside [left, right] ────────────────────*/
    std::vector<Band> find_bands_in_region(size_t left_bound, size_t right_bound) {
        std::vector<Band> bands;

        const size_t n = pairings_.size();
        if (right_bound >= n) {
            right_bound = n - 1;
        }

        for (size_t i = left_bound; i <= right_bound; ++i) {
            /* skip anything that is   – already used
             *                         – unpaired
             *                         – closing half
             *                         – pairs outside this region */
            if (done_[i] || pairings_[i] == NULL_INDEX || pairings_[i] < i ||
                pairings_[i] > right_bound)
                continue;

            size_t j = pairings_[i];
            size_t i_prime = i;
            size_t j_prime = j;

            // walks the stem until last base pair in the given region (finds i` and j`)
            while (extend_stem(i_prime, j_prime, left_bound, right_bound)) {
            }

            // stores entire band
            bands.push_back(Band{i, i_prime, j_prime, j, pairings_});

            /* mark every position that belongs to this band */
            for (size_t k = i; k <= i_prime; ++k) {
                done_[k] = true;
                if (pairings_[k] != NULL_INDEX) done_[pairings_[k]] = true;
            }
            for (size_t k = j_prime; k <= j; ++k) {
                done_[k] = true;
                if (pairings_[k] != NULL_INDEX) done_[pairings_[k]] = true;
            }

            i = i_prime;  // fast-forward
        }
        return bands;
    }

    std::vector<Band> find_bands_in_region(size_t start) {
        return find_bands_in_region(start, pairings_.size() - 1);
    }

    std::vector<Band> find_bands_in_region() {
        return find_bands_in_region(0, pairings_.size() - 1);
    }
};

}  // namespace knotergy