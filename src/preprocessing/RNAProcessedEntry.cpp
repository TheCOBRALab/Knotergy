#include "RNAProcessedEntry.hpp"

namespace knotergy {

RNAProcessedEntry::RNAProcessedEntry(const RNAEntry& rna) : rna_{rna} {
    if (rna_.sequence.size() != rna_.structure.size()) {
        THROW_ERROR("RNA Sequence & Structure's lengths are miss-matched. Sequence length: " +
                    std::to_string(rna_.sequence.size()) +
                    " Structure: " + std::to_string(rna_.structure.size()));
    }
    pairings_ = compute_pairings();
    closed_regions_ = compute_closed_regions();
    closed_region_pairings_ = compute_closed_region_pairings();
    unpaired_prefix_sum_ = compute_unpaired_counts();
}

RNAProcessedEntry::RNAProcessedEntry(std::string name, std::string sequence,
                                     std::string structure) {
    RNAProcessedEntry(RNAEntry{name, sequence, structure});
}

RNAProcessedEntry::RNAProcessedEntry(std::string sequence, std::string structure) {
    RNAProcessedEntry(RNAEntry{sequence, structure});
}

const std::string& RNAProcessedEntry::get_name() const {
    return rna_.name;
}

const std::string& RNAProcessedEntry::get_sequence() const {
    return rna_.sequence;
}

const std::string& RNAProcessedEntry::get_structure() const {
    return rna_.structure;
}

const std::vector<size_t>& RNAProcessedEntry::get_pairings() const {
    return pairings_;
}

const std::vector<size_t>& RNAProcessedEntry::get_closed_region_pairings() const {
    return closed_region_pairings_;
}

const std::vector<ClosedRegion>& RNAProcessedEntry::get_closed_regions() const {
    return closed_regions_;
}

size_t RNAProcessedEntry::size() const {
    return rna_.size();
}

int RNAProcessedEntry::get_unpaired_count(size_t from, size_t to) const {
    if (from >= unpaired_prefix_sum_.size() || to > unpaired_prefix_sum_.size()) {
        throw std::out_of_range("Index out of range in get_unpaired_count");
    }

    if (from >= to) return 0;

    // Return the difference in unpaired counts between the two indices
    return unpaired_prefix_sum_[to] - unpaired_prefix_sum_[from];
}

int RNAProcessedEntry::get_unpaired_count(ClosedRegion cr) const {
    return get_unpaired_count(cr.begin, cr.end);
}

std::vector<size_t> RNAProcessedEntry::compute_pairings() {
    std::unordered_map<char, char> open_to_close = {{'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}};

    // close to open is the opposite of open_to_close
    // e.g. {'(', ')'} -> {')', '('}
    std::unordered_map<char, char> close_to_open;
    close_to_open.reserve(open_to_close.size());
    for (const auto& pair : open_to_close) {
        close_to_open.emplace(pair.second, pair.first);
    }

    std::unordered_map<char, char> valid_pairings = {{'A', 'U'}, {'U', 'A'}, {'G', 'C'},
                                                     {'C', 'G'}, {'G', 'U'}, {'U', 'G'}};

    std::vector<size_t> pairings(rna_.size(), NULL_INDEX);
    std::unordered_map<char, std::vector<size_t>> stacks;
    for (size_t i = 0; i < rna_.size(); i++) {
        char c = rna_.structure[i];

        if (c == '.') continue;

        // if open bracket
        if (open_to_close.count(c)) {
            stacks[c].push_back(i);
            continue;
        }

        // if closing bracket
        if (close_to_open.count(c)) {
            char open = close_to_open[c];
            std::vector<size_t>& stack = stacks[open];

            if (stacks[open].empty()) {
                THROW_ERROR("Invalid RNA structure: Bracket: '" + std::string(1, c) +
                            "' at index: " + std::to_string(i) + " was never opened");
            }
            size_t j = stack.back();
            stack.pop_back();
            pairings[i] = j;
            pairings[j] = i;
            if (valid_pairings[rna_.sequence[i]] != rna_.sequence[j]) {
                std::cerr << "Warning: Base Pair '" + std::string(1, rna_.sequence[i]) +
                                 "' can't pair with " + std::string(1, rna_.sequence[j]) + '\''
                          << std::endl;
            }
            continue;
        }

        THROW_ERROR("Invalid RNA structure: Invalid character '" + std::string(1, c) +
                    "' in RNA structure (position " + std::to_string(i) + ')');
    }

    for (auto& [open, stack] : stacks) {
        if (!stack.empty()) {
            THROW_ERROR("Invalid RNA structure: opening bracket '" + std::string(1, open) +
                        "' at index: " + std::to_string(stack.back()) + " was not closed");
        }
    }

    return pairings;
}

std::vector<ClosedRegion> RNAProcessedEntry::compute_closed_regions() {
    std::vector<ClosedRegion> closed_regions;
    std::stack<ClosedRegion> stack;
    const size_t n = rna_.size();

    for (size_t i = 0; i < n; ++i) {
        size_t bp = pairings_[i];
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
std::vector<size_t> RNAProcessedEntry::compute_closed_region_pairings() {
    std::vector<size_t> closed_region_pairings(rna_.size(), NULL_INDEX);
    for (ClosedRegion cr : closed_regions_){
        closed_region_pairings[cr.begin] = cr.end;
        closed_region_pairings[cr.end] = cr.begin;
    }
    return closed_region_pairings;
}

std::vector<int> RNAProcessedEntry::compute_unpaired_counts() {
    int count = 0;
    size_t n = rna_.size();
    std::vector<int> unpaired_count_list;

    unpaired_count_list.assign(n + 1, 0);
    for (size_t i = 0; i < n; ++i) {
        count += (pairings_[i] == NULL_INDEX);
        unpaired_count_list[i + 1] = count;
    }
    return unpaired_count_list;
};

}  // namespace knotergy
