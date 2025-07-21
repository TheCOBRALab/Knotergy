#include "RNAProcessedEntry.hpp"

namespace knotergy {

RNAProcessedEntry::RNAProcessedEntry(const RNAEntry& rna) : rna_{rna} {
    pairings_ = compute_pairings();
    closed_regions_ = compute_closed_regions();
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
    std::stack<size_t> brackets;
    std::stack<size_t> pseudoknots;
    size_t j;
    std::vector<size_t> pairings(rna_.size(), NULL_INDEX);

    for (size_t i = 0; i < rna_.size(); i++) {
        switch (rna_.structure[i]) {
            case '.':
                break;
            case '(':
                brackets.push(i);
                break;
            case '[':
                pseudoknots.push(i);
                break;
            case ')':
                if (brackets.empty()) {
                    throw std::runtime_error("Structure in RNAEntry is invalid. \nSequence: " +
                                             rna_.sequence + "\nStructure: " + rna_.structure);
                }
                j = brackets.top();
                brackets.pop();
                pairings[i] = j;
                pairings[j] = i;
                break;
            case ']':
                if (pseudoknots.empty()) {
                    throw std::runtime_error("Structure in RNAEntry is invalid. \nSequence: " +
                                             rna_.sequence + "\nStructure: " + rna_.structure);
                }
                j = pseudoknots.top();
                pseudoknots.pop();
                pairings[i] = j;
                pairings[j] = i;
                break;
            default:
                throw std::runtime_error(
                    "Character in RNAEntry's structure is invalid. \nInvalid Character: " +
                    std::string(1, rna_.structure[i]) + "\nSequence: " + rna_.sequence +
                    "\nStructure: " + rna_.structure);
        }
    }
    if (!brackets.empty() || !pseudoknots.empty()) {
        throw std::runtime_error("Unmatched opening brackets in RNA structure.\nSequence: " +
                                 rna_.sequence + "\nStructure: " + rna_.structure);
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
