#include "BandBuilder.hpp"

#include "utils/common.hpp"

namespace knotergy {

Band BandBuilder::construct_band(size_t lb, size_t li, size_t ri, size_t rb,
                                 const std::vector<size_t>& pair_table,
                                 const std::vector<size_t>& cr_pair_table) {
    // ------------- Validate band structure and pairing -------------
    if (lb >= pair_table.size() || li >= pair_table.size() || ri >= pair_table.size() ||
        rb >= pair_table.size()) {
        THROW_ERROR("Band boundary indices must be within the bounds of the sequence.\n");
    }

    if (li < lb || ri <= li || rb < ri) {
        THROW_ERROR(
            "Band boundaries must satisfy: left border < left inner < right inner < right "
            "border.\n");
    }

    if (pair_table[lb] != rb || pair_table[li] != ri) {
        THROW_ERROR("Band boundary indices must be correctly paired in the pair table.\n");
    };

    // ------------- Build and populate base pairs -------------
    int                   child_count = 0;
    std::vector<BasePair> base_pairs =
        find_base_pairs_left_scan(lb, li, ri, rb, pair_table, cr_pair_table, child_count);
    populate_right_arm_children(base_pairs, ri, rb, cr_pair_table, child_count);

    return Band(lb, li, ri, rb, std::move(base_pairs), child_count);
}

Band BandBuilder::construct_band(BandBounds bounds, const std::vector<size_t>& pair_table,
                                 const std::vector<size_t>& cr_pair_table) {
    return construct_band(bounds.left_border, bounds.left_inner, bounds.right_inner,
                          bounds.right_border, pair_table, cr_pair_table);
}

std::vector<BasePair> BandBuilder::find_base_pairs_left_scan(
    size_t lb, size_t li, size_t ri, size_t rb, const std::vector<size_t>& pair_table,
    const std::vector<size_t>& cr_pair_table, int& child_count) {
    std::vector<BasePair> base_pairs;
    base_pairs.reserve(std::min(rb - ri, li - lb) + 1);  // Max possible base pairs in the band

    base_pairs.emplace_back(lb, pair_table[lb]);

    for (size_t idx = lb + 1; idx <= li; ++idx) {
        // Skip closed region and add it as a child of the current base pair
        if (cr_pair_table[idx] != NULL_INDEX) {
            base_pairs.back().children.emplace_back(idx, cr_pair_table[idx]);
            idx = cr_pair_table[idx];
            ++child_count;
            continue;
        }

        // Check if the current index is a base pair that belongs to the band
        size_t paired = pair_table[idx];
        if (paired >= ri && paired <= rb) {
            base_pairs.emplace_back(idx, paired);
        }
    }

    return base_pairs;
}

void BandBuilder::populate_right_arm_children(std::vector<BasePair>& base_pairs, size_t ri,
                                              size_t rb, const std::vector<size_t>& cr_pair_table,
                                              int& child_count) {
    if (base_pairs.empty()) {
        return;
    }

    size_t current_bp_idx = 0;  // Used to track which base pair we are currently adding children to

    size_t next_bp_right_border = base_pairs.size() > 1 ? base_pairs[1].j : NULL_INDEX;

    // Scans from right border towards right inner.
    // next_bp_right_border tracks the right base of the next base pair in the band
    for (size_t idx = rb - 1; idx > ri; --idx) {
        BasePair& current_bp = base_pairs[current_bp_idx];

        // Keeps track of the current base pair we are adding children to.
        if (idx == next_bp_right_border) {
            ++current_bp_idx;

            if (base_pairs.size() > current_bp_idx + 1) {
                next_bp_right_border = base_pairs[current_bp_idx + 1].j;
            } else {
                next_bp_right_border = NULL_INDEX;
            }

            continue;
        }

        // Adds children to current base pair
        if (cr_pair_table[idx] != NULL_INDEX) {
            size_t right = idx;
            size_t left  = cr_pair_table[idx];

            current_bp.children.emplace_back(left, right);
            idx = left;

            ++child_count;

            if (idx == 0) {
                break;
            }

            continue;
        }
    }
};

}  // namespace knotergy