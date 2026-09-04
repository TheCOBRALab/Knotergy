#include "BandBuilder.hpp"

#include "utils/common.hpp"

namespace knotergy {

Band BandBuilder::construct_band(std::size_t lb, std::size_t li, std::size_t ri, std::size_t rb,
                                 const std::vector<std::size_t>& pair_table,
                                 const std::vector<LoopNode*>& node_table) {
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

    if (node_table.size() != pair_table.size()) {
        THROW_ERROR("Loop-node lookup size must match pair table size.");
    }

    // ------------- Build and populate base pairs -------------
    int child_count = 0;
    std::vector<PKBasePair> base_pairs =
        find_base_pairs_left_scan(lb, li, ri, rb, pair_table, node_table, child_count);
    populate_right_arm_children(base_pairs, ri, rb, node_table, child_count);

    return Band(lb, li, ri, rb, std::move(base_pairs), child_count);
}

Band BandBuilder::construct_band(BandBounds bounds, const std::vector<std::size_t>& pair_table,
                                 const std::vector<LoopNode*>& node_table) {
    return construct_band(bounds.left_border, bounds.left_inner, bounds.right_inner,
                          bounds.right_border, pair_table, node_table);
}

std::vector<PKBasePair> BandBuilder::find_base_pairs_left_scan(
    std::size_t lb, std::size_t li, std::size_t ri, std::size_t rb,
    const std::vector<std::size_t>& pair_table, const std::vector<LoopNode*>& node_table,
    int& child_count) {
    std::vector<PKBasePair> base_pairs;
    base_pairs.reserve(std::min(rb - ri, li - lb) + 1);  // Max possible base pairs in the band

    base_pairs.emplace_back(lb, pair_table[lb]);

    for (std::size_t idx = lb + 1; idx <= li; ++idx) {
        // Skip closed region and add it as a child of the current base pair
        if (node_table[idx] != nullptr) {
            node_table[idx]->pseudo_type = PseudoNestedType::WithinBand;
            base_pairs.back().children.emplace_back(idx, node_table[idx]->end);
            idx = node_table[idx]->end;
            ++child_count;
            continue;
        }

        // Check if the current index is a base pair that belongs to the band
        std::size_t paired = pair_table[idx];
        if (paired >= ri && paired <= rb) {
            base_pairs.emplace_back(idx, paired);
        }
    }

    return base_pairs;
}

void BandBuilder::populate_right_arm_children(std::vector<PKBasePair>& base_pairs, std::size_t ri,
                                              std::size_t rb,
                                              const std::vector<LoopNode*>& node_table,
                                              int& child_count) {
    if (base_pairs.empty()) {
        return;
    }

    std::size_t current_bp_idx =
        0;  // Used to track which base pair we are currently adding children to

    std::size_t next_bp_right_border = base_pairs.size() > 1 ? base_pairs[1].j : NULL_INDEX;

    // Scans from right border towards right inner.
    // next_bp_right_border tracks the right base of the next base pair in the band
    for (std::size_t idx = rb - 1; idx > ri; --idx) {
        PKBasePair& current_bp = base_pairs[current_bp_idx];

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
        if (node_table[idx] != nullptr) {
            std::size_t right = idx;
            std::size_t left = node_table[idx]->begin;

            node_table[idx]->pseudo_type = PseudoNestedType::WithinBand;

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