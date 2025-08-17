#include "RNAProcessor.hpp"

namespace knotergy {
ProcessedRNAEntry RNAProcessor::process_rna(RNAEntry rna) {
    std::vector<size_t> pairings = compute_pairings(rna);
    std::vector<ClosedRegion> closed_regions = compute_closed_regions(pairings);
    std::vector<size_t> cr_pairings = compute_closed_regions_pairings(closed_regions, rna.size());
    std::vector<int> unpaired_prefix_sum = compute_unpaired_counts(pairings);

    return ProcessedRNAEntry{std::move(rna), std::move(pairings), std::move(closed_regions),
                             std::move(cr_pairings), std::move(unpaired_prefix_sum)};
};

std::vector<size_t> RNAProcessor::compute_pairings(RNAEntry& rna) {
    std::unordered_map<char, char> open_to_close = {{'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}};
    std::unordered_map<char, char> valid_pairings = {{'A', 'U'}, {'U', 'A'}, {'G', 'C'},
                                                     {'C', 'G'}, {'G', 'U'}, {'U', 'G'}};

    // close to open is the opposite of open_to_close
    // e.g. {'(', ')'} -> {')', '('}
    std::unordered_map<char, char> close_to_open;
    close_to_open.reserve(open_to_close.size());
    for (const auto& pair : open_to_close) {
        close_to_open.emplace(pair.second, pair.first);
    }

    // pre-allocate pairings
    std::vector<size_t> pairings(rna.size(), NULL_INDEX);

    // 1 stack for each open/close pair
    std::unordered_map<char, std::vector<size_t>> stacks;

    for (size_t i = 0; i < rna.size(); i++) {
        char c = rna.structure[i];

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
            if (valid_pairings[rna.sequence[i]] != rna.sequence[j]) {
                std::cerr << "Warning: Base Pair '" + std::string(1, rna.sequence[i]) +
                                 "' can't pair with " + std::string(1, rna.sequence[j]) + '\''
                          << std::endl;
            }
            continue;
        }

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

std::vector<ClosedRegion> RNAProcessor::compute_closed_regions(std::vector<size_t>& pairings) {
    std::vector<ClosedRegion> closed_regions;
    std::stack<ClosedRegion> stack;
    const size_t n = pairings.size();

    for (size_t i = 0; i < n; ++i) {
        size_t bp = pairings[i];
        if (bp == NULL_INDEX) continue;  // unpaired

        // ───── OPENING BASE: i < bp ────────────────────────────
        if (i < bp) {
            stack.push({i, bp});
            continue;
        }

        // ───── CLOSING BASE: bp < i ────────────────────────────
        size_t largest_right = i;  // rightmost boundary seen

        // if crossing (pseudoknotted), find right end of closed region
        while (!stack.empty() && stack.top().begin > bp) {
            largest_right = std::max(largest_right, stack.top().end);
            stack.pop();
        }

        if (stack.empty()) continue;  // if unbalanced (should never happen)

        // extend region if needed
        stack.top().end = std::max(largest_right, stack.top().end);

        // region finished?
        if (i == stack.top().end) {
            closed_regions.push_back(stack.top());
            stack.pop();
        }
    }
    return closed_regions;
}

// ([...)] = 6, -1, -1, -1, -1, -1, 0
std::vector<size_t> RNAProcessor::compute_closed_regions_pairings(
    std::vector<ClosedRegion>& closed_regions, const size_t& rna_size) {
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

}  // namespace knotergy