#include "RNAProcessor.hpp"

#include <unordered_map>
#include <unordered_set>
namespace knotergy {
ProcessedRNAEntry RNAProcessor::process_rna(
    RNAEntry rna, const std::vector<modified_base_param>& modified_params) {
    std::vector<std::string_view> mod_sequence = ProcessedRNAEntry::compute_modified_sequence_views(rna.sequence);

    // Ensure Sequence & Structure are the same length
    if (mod_sequence.size() != rna.structure.size()) {
        THROW_ERROR("Modified sequence length does not match RNA length\nSequence length: " +
                    std::to_string(mod_sequence.size()) +
                    "\nRNA length: " + std::to_string(rna.structure.size()));
    }

    bool has_modified_bases = false;
    std::string unmodified_sequence = compute_unmodified_sequence(mod_sequence, modified_params, rna.size(), has_modified_bases);
    std::vector<size_t> pairings = compute_pairings(rna, unmodified_sequence, mod_sequence);
    std::vector<ClosedRegion> closed_regions = compute_closed_regions(pairings);
    std::vector<size_t> cr_pairings = compute_cr_pairings(closed_regions, rna.size());
    std::vector<int> unpaired_prefix_sum = compute_unpaired_counts(pairings);

    return ProcessedRNAEntry{std::move(rna),
                             std::move(unmodified_sequence),
                             std::move(pairings),
                             std::move(closed_regions),
                             std::move(cr_pairings),
                             std::move(unpaired_prefix_sum),
                             has_modified_bases};
};

std::vector<size_t> RNAProcessor::compute_pairings(
    const std::string& structure, const std::string& unmodified_sequence,
    const std::vector<std::string_view>& mod_sequence) {
    const std::unordered_map<char, char> open_to_close = {
        {'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}, {'A', 'a'}, {'B', 'b'},
        {'C', 'c'}, {'D', 'd'}, {'E', 'e'}, {'F', 'f'}, {'G', 'g'}, {'H', 'h'},
        {'I', 'i'}, {'J', 'j'}, {'K', 'k'}, {'L', 'l'}, {'M', 'm'}, {'N', 'n'},
        {'O', 'o'}, {'P', 'p'}, {'Q', 'q'}, {'R', 'r'}, {'S', 's'}, {'T', 't'},
        {'U', 'u'}, {'V', 'v'}, {'W', 'w'}, {'X', 'x'}, {'Y', 'y'}, {'Z', 'z'}};

    const std::unordered_map<char, std::unordered_set<char>> valid_pairings = {
        {'A', {'U'}}, {'U', {'A', 'G'}},      {'G', {'C', 'U'}},
        {'C', {'G'}}, {'T', {'A', 'T', 'G'}}, {'N', {}}};

    // close_to_open is the opposite of open_to_close
    // e.g. {'(', ')'} -> {')', '('}
    std::unordered_map<char, char> close_to_open;
    close_to_open.reserve(open_to_close.size());
    for (const auto& pair : open_to_close) {
        close_to_open.emplace(pair.second, pair.first);
    }

    // pre-allocate pairings
    std::vector<size_t> pairings(structure.size(), NULL_INDEX);

    // 1 stack for each open/close pair
    std::unordered_map<char, std::vector<size_t>> stacks;

    for (size_t i = 0; i < structure.size(); i++) {
        char c = structure[i];

        // un-paired
        if (c == '.') continue;

        // if open bracket, push to it's stack
        if (open_to_close.count(c)) {
            stacks[c].push_back(i);
            continue;
        }

        // if closing bracket
        if (close_to_open.count(c)) {
            char open = close_to_open[c];
            std::vector<size_t>& stack = stacks[open];

            // if closed without being opened
            if (stacks[open].empty()) {
                THROW_ERROR("Invalid RNA structure: Bracket: '" + std::string(1, c) +
                            "' at index: " + std::to_string(i) + " was never opened");
            }

            // pair top of stack with current basepair
            size_t j = stack.back();
            stack.pop_back();
            pairings[i] = j;
            pairings[j] = i;

            // check if they're a valid pair
            if (unmodified_sequence.size() == structure.size() && !valid_pairings.at(unmodified_sequence[j]).count(unmodified_sequence[i])) {
                // if modified sequence is provided, use that for warning
                if (mod_sequence.size() == structure.size()) { 

                    std::cerr << "Warning: Base Pair '" + std::string(mod_sequence[i]) +
                                     "' can't pair with '" + std::string(mod_sequence[j]) +
                                     "' at indices " + std::to_string(j) + ", " + std::to_string(i)
                              << std::endl;
                } else {  // If modified sequence is not provided, fall back to unmodified sequence
                    std::cerr << "Warning: Base Pair '" + std::string(1, unmodified_sequence[i]) +
                                 "' can't pair with '" + std::string(1, unmodified_sequence[j]) +
                                 "' at indices " + std::to_string(j) + ", " + std::to_string(i)
                                 << std::endl;
                }
            }
            continue;
        }

        // Not open, not close, not unpaired -> invalid character
        THROW_ERROR("Invalid RNA structure: Invalid character '" + std::string(1, c) +
                    "' in RNA structure (position " + std::to_string(i) + ')');
    }

    // check if all open brackets were closed
    for (auto& [open, stack] : stacks) {
        if (!stack.empty()) {
            THROW_ERROR("Invalid RNA structure: opening bracket '" + std::string(1, open) +
                        "' at index: " + std::to_string(stack.back()) + " was not closed");
        }
    }

    return pairings;
}

std::vector<size_t> RNAProcessor::compute_pairings(
    const RNAEntry& rna, const std::string& unmodified_sequence,
    const std::vector<std::string_view>& mod_sequence) {
    return compute_pairings(rna.structure, unmodified_sequence, mod_sequence);
}

// Computes all closed regions from a pairing list using an interval-merging stack algorithm.
// Assumes properly balanced pairings.
std::vector<ClosedRegion> RNAProcessor::compute_closed_regions(
    const std::vector<size_t>& pairings) {
    std::vector<ClosedRegion> closed_regions;
    std::stack<ClosedRegion> stack;
    const size_t n = pairings.size();

    for (size_t i = 0; i < n; ++i) {
        size_t paired_idx = pairings[i];
        if (paired_idx == NULL_INDEX) continue;

        if (i < paired_idx) {
            // Opening base: push raw pair; its right boundary may be extended later.
            stack.emplace(i, paired_idx);
            continue;
        }

        const size_t open_idx = paired_idx;
        size_t largest_right = i;

        // Merge any nested regions whose start lies within (open, i]
        while (!stack.empty() && (stack.top().begin > open_idx)) {
            largest_right = std::max(largest_right, stack.top().end);
            stack.pop();
        }

        if (stack.empty()) THROW_ERROR("Unbalanced pairings\n");

        // Extend region if any nested intervals were merged.
        stack.top().end = std::max(largest_right, stack.top().end);

        // If this closing base ends the region, move it to result.
        if (i == stack.top().end) {
            closed_regions.push_back(stack.top());
            stack.pop();
        }
    }
    return closed_regions;
}

// ([...)] = 6, -1, -1, -1, -1, -1, 0
std::vector<size_t> RNAProcessor::compute_cr_pairings(
    const std::vector<ClosedRegion>& closed_regions, const size_t& rna_size) {
    std::vector<size_t> closed_regions_pairings(rna_size, NULL_INDEX);
    for (ClosedRegion cr : closed_regions) {
        if (cr.end >= rna_size) THROW_ERROR("rna_size is too small");
        closed_regions_pairings[cr.begin] = cr.end;
        closed_regions_pairings[cr.end] = cr.begin;
    }
    return closed_regions_pairings;
}

std::vector<int> RNAProcessor::compute_unpaired_counts(const std::vector<size_t>& pairings) {
    int count = 0;
    size_t n = pairings.size();
    std::vector<int> unpaired_prefix_sum;

    unpaired_prefix_sum.assign(n + 1, 0);
    for (size_t i = 0; i < n; ++i) {
        count += (pairings[i] == NULL_INDEX);
        unpaired_prefix_sum[i + 1] = count;
    }
    return unpaired_prefix_sum;
};

std::string RNAProcessor::compute_unmodified_sequence(
    const std::vector<std::string_view>& modified_sequence_views,
    const std::vector<modified_base_param>& params, const size_t rna_length,
    bool& has_modified_bases) {
    std::string unmodified_sequence;
    unmodified_sequence.reserve(rna_length);

    // map modified base to unmodified base for quick lookup
    std::unordered_map<std::string_view, std::string> mod_to_unmod;
    mod_to_unmod.reserve(params.size());
    for (const modified_base_param& param : params) {
        if (param.unmodified_base.empty()) {
            THROW_ERROR("Modified base '" + param.modified_base +
                        "' has empty unmodified mapping.");
        }

        mod_to_unmod.emplace(param.modified_base, param.unmodified_base);
    }

    // Convert modified sequence to unmodified sequence
    for (const std::string_view& mod_base : modified_sequence_views) {
        if (mod_base.empty()) {
            THROW_ERROR("Empty base in modified sequence");
        }

        if (is_unmod_base(mod_base)) {
            unmodified_sequence.append(mod_base);
            continue;
        } else {
            has_modified_bases = true;
        }

        // lookup modified -> unmodified
        auto it = mod_to_unmod.find(mod_base);
        if (it != mod_to_unmod.end()) {
            unmodified_sequence += it->second;
        } else {
            THROW_ERROR("Base '" + std::string(mod_base) + "' at index " +
                        std::to_string(unmodified_sequence.size()) + " is not a valid base (not found in modified params and is not a standard unmodified base)");
        }
    }

    return unmodified_sequence;
}

// Check if a base is unmodified (using lookup table for speed)
#include <iostream>

bool RNAProcessor::is_unmod_base(const std::string_view& b) {
    if (b.size() != 1) return false; // Multi-character "base" are considered modified (e.g. ❤️‍🩹)
    return unmod_lookup[(unsigned char)b[0]] != 0;
}


}  // namespace knotergy