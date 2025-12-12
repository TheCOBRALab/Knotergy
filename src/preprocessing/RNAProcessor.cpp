#include "RNAProcessor.hpp"
#include <unordered_map>
#include <unordered_set>

namespace knotergy {
ProcessedRNAEntry RNAProcessor::process_rna(RNAEntry rna) {
    std::vector<std::string_view> mod_sequence = compute_modified_sequence_views(rna);
    std::vector<size_t> pairings = compute_pairings(rna);
    std::vector<ClosedRegion> closed_regions = compute_closed_regions(pairings);
    std::vector<size_t> cr_pairings = compute_cr_pairings(closed_regions, rna.size());
    std::vector<int> unpaired_prefix_sum = compute_unpaired_counts(pairings);

    return ProcessedRNAEntry{std::move(rna), std::move(mod_sequence), std::move(pairings),
                             std::move(closed_regions), std::move(cr_pairings),
                             std::move(unpaired_prefix_sum)};
};

std::vector<size_t> RNAProcessor::compute_pairings(const RNAEntry& rna) {
    std::unordered_map<char, char> open_to_close = {
        {'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'},
        {'A', 'a'}, {'B', 'b'}, {'C', 'c'}, {'D', 'd'}
    };
    std::unordered_map<char, std::unordered_set<char>> valid_pairings = {
        {'A', {'U'}},
        {'U', {'A','G'}},
        {'G', {'C','U'}},
        {'C', {'G'}}
    };

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
            if (!valid_pairings[rna.sequence[j]].count(rna.sequence[i])) {
                std::cerr << "Warning: Base Pair '" + std::string(1, rna.sequence[i]) +
                                 "' can't pair with '" + std::string(1, rna.sequence[j]) + "' at indices " + std::to_string(j) + ", " + std::to_string(i) << std::endl;
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

// Computes all closed regions from a pairing list using an interval-merging stack algorithm.
// Assumes properly balanced pairings.
std::vector<ClosedRegion> RNAProcessor::compute_closed_regions(const std::vector<size_t>& pairings) {
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

std::vector<std::string_view> RNAProcessor::compute_modified_sequence_views(const RNAEntry& rna) {
    std::vector<std::string_view> mod_sequence;
    mod_sequence.reserve(rna.sequence.size());
    for (size_t head = 0; head < rna.sequence.size();) {
        const size_t len = bytes_to_read(rna.sequence, head);
        if (len == 0) {
            ++head;
            continue;
        }
        mod_sequence.emplace_back(std::string_view(&rna.sequence[head], len));
        head += len;
    }

    // for (const std::string_view& sv : mod_sequence) {
    //     std::cout << sv << '\n';
    // }
    return mod_sequence;
};


// Checks if the substring of s starting at pos matches the given byte pattern
bool RNAProcessor::check_utf8_pattern(const std::string& s,
                                 size_t pos,
                                 const std::vector<unsigned char>& pattern) noexcept
{
    const size_t patternLen = pattern.size();
    if (pos + patternLen > s.size()) return false;
    for (size_t i = 0; i < patternLen; ++i)
        if (static_cast<unsigned char>(s[pos + i]) != pattern[i]) return false;
    return true;
}

// Returns number of bytes to read for the UTF-8 character starting at sequence[index]
size_t RNAProcessor::bytes_to_read(const std::string& sequence, size_t index) noexcept {
    const size_t n = sequence.size();
    if (index >= n) return 0;
    const unsigned char first = static_cast<unsigned char>(sequence[index]);

    // Determine the number of bytes to read based on the first byte
    size_t bytesToRead =
        ((first & 0x80u) == 0x00u) ? 1 :
        ((first & 0xE0u) == 0xC0u) ? 2 :
        ((first & 0xF0u) == 0xE0u) ? 3 :
        ((first & 0xF8u) == 0xF0u) ? 4 : 1;

    // Ensure we don't read past the end of the string
    if (bytesToRead > n - index) return n - index;

    size_t pos = index + bytesToRead;

    // Check for variation selectors and zero-width joiners
    static const std::vector<unsigned char> ZWJ_BYTES  = {0xE2u, 0x80u, 0x8Du}; // U+200D
    static const std::vector<unsigned char> VS15_BYTES = {0xEFu, 0xB8u, 0x8Eu}; // U+FE0E
    static const std::vector<unsigned char> VS16_BYTES = {0xEFu, 0xB8u, 0x8Fu}; // U+FE0F
    
    // Skin tones e.g. 👍🏿 = 👍 + 🏿
    if (index + bytesToRead + 4 <= n &&
        static_cast<unsigned char>(sequence[index + bytesToRead]    )     == 0xF0u && 
        static_cast<unsigned char>(sequence[index + bytesToRead + 1])     == 0x9Fu && 
        static_cast<unsigned char>(sequence[index + bytesToRead + 2])     == 0x8Fu &&
        static_cast<unsigned char>(sequence[index + bytesToRead + 3])     >= 0xBBu &&
        static_cast<unsigned char>(sequence[index + bytesToRead + 3])     <= 0xBFu)
    {
        bytesToRead += 4;
        pos += 4;
    }
    
    // Variation Selectors (VS15, VS16) e.g. 😐︎ = 😐 + VS15; ❤️ = ❤ + VS16
    if (check_utf8_pattern(sequence, pos, VS16_BYTES) ||
        check_utf8_pattern(sequence, pos, VS15_BYTES)){
        bytesToRead += VS16_BYTES.size();
        pos += VS16_BYTES.size();
    }

    // Zero-width joiner (ZWJ) e.g. ❤️‍🩹 = ❤️ + ZWJ + 🩹
    if (check_utf8_pattern(sequence, pos, ZWJ_BYTES)) {
        bytesToRead += ZWJ_BYTES.size();
        bytesToRead += bytes_to_read(sequence, pos + ZWJ_BYTES.size());
        return bytesToRead;
    }

    if (bytesToRead > n - index) return n - index;
    return bytesToRead;
}

}  // namespace knotergy