#include "RNAProcessor.hpp"

#include <unordered_map>
#include <unordered_set>
namespace knotergy {
ProcessedRNAEntry RNAProcessor::process_rna(RNAEntry rna, const all_mod_params& modified_params) {
    bool has_modified_bases = false;
    std::string unmodified_sequence;
    std::vector<size_t> pairings;

    if (!modified_params.empty()) {
        std::vector<std::string_view> mod_sequence =
            ProcessedRNAEntry::compute_modified_sequence_views(rna.sequence, rna.structure);

        if (mod_sequence.size() != rna.structure.size()) {
            THROW_ERROR("Sequence length does not match RNA length\nSequence length: " +
                        std::to_string(mod_sequence.size()) +
                        "\nRNA length: " + std::to_string(rna.structure.size()));
        }

        unmodified_sequence = compute_unmodified_sequence(mod_sequence, modified_params, rna.size(),
                                                          has_modified_bases);
        pairings = compute_pairings(rna, unmodified_sequence, mod_sequence);
    }
    // If no modified base parameters provided, skip modified sequence processing
    else {
        // Validate that the sequence only contains unmodified bases, and warn if it contains
        // modified bases but no parameters are provided
        for (char base : rna.sequence) {
            if (unmod_lookup[(unsigned char) base] == 0) {
                THROW_ERROR("RNA sequence contains invalid base: '" + std::string(1, base) + "'");
                break;
            }
        }

        unmodified_sequence = rna.sequence;
        pairings = compute_pairings(rna, unmodified_sequence);
    }

    std::vector<ClosedRegion> closed_regions = compute_closed_regions(pairings);
    std::vector<size_t> cr_pairings = compute_cr_pairings(closed_regions, rna.size());
    std::vector<int> unpaired_prefix_sum = compute_unpaired_counts(pairings);

    return ProcessedRNAEntry{std::move(rna),         std::move(unmodified_sequence),
                             std::move(pairings),    std::move(closed_regions),
                             std::move(cr_pairings), std::move(unpaired_prefix_sum),
                             has_modified_bases};
};

std::vector<size_t> RNAProcessor::compute_pairings(
    const std::string& structure, const std::string& unmodified_sequence,
    const std::vector<std::string_view>& mod_sequence) {
    constexpr unsigned char MAX_CHAR = 128;
    constexpr unsigned char INVALID_CHAR = '\0';

    // Lookup tables for bracket matching. Indexed by ASCII character code.
    std::array<char, MAX_CHAR> open_to_close{};  // open('(') -> close(')')
    std::array<char, MAX_CHAR> close_to_open{};  // close(')') -> open('(')

    // Adds valid bracket pairs to the lookup tables
    auto add_bracket = [&](char open, char close) {
        open_to_close[static_cast<unsigned char>(open)] = close;
        close_to_open[static_cast<unsigned char>(close)] = open;
    };

    add_bracket('(', ')');
    add_bracket('[', ']');
    add_bracket('{', '}');
    add_bracket('<', '>');

    // Adds A-Z and a-z as valid brackets. e.g. open('A') -> close('a')
    for (char c = 'A'; c <= 'Z'; ++c) {
        add_bracket(c, static_cast<char>(std::tolower(c)));
    }

    // Function to check if two bases can pair according to RNA base-pairing rules.
    auto can_pair = [](char left, char right) {
        switch (left) {
            case 'A':
                return right == 'U';
            case 'U':
                return right == 'A' || right == 'G';
            case 'G':
                return right == 'C' || right == 'U';
            case 'C':
                return right == 'G';
            case 'T':
                return right == 'A' || right == 'G';
            case 'N':
                return false;  // Should not pair with anything
            default:
                return false;
        }
    };

    const size_t n = structure.size();

    // sequence and mod_sequence are optional inputs
    // If provided, it will check for valid base pairs, and warn for invalid pairs.
    const bool check_pairs = unmodified_sequence.size() == n;
    const bool have_mod_seq = mod_sequence.size() == n;

    // Initialize all pairings to NULL_INDEX (unpaired).
    std::vector<size_t> pairings(n, NULL_INDEX);

    // Stacks for each type of opening bracket. Indexed by ASCII character code.
    std::array<std::vector<size_t>, MAX_CHAR> stacks;

    for (size_t i = 0; i < n; ++i) {
        const unsigned char c = static_cast<unsigned char>(structure[i]);

        // Skip unpaired bases
        if (c == '.') {
            continue;
        }
        // If it's an opening bracket, push its index onto the corresponding stack
        if (c < MAX_CHAR && open_to_close[c] != INVALID_CHAR) {
            stacks[c].push_back(i);
            continue;
        }

        // If it's a closing bracket, pop from the corresponding opening stack and record the pair
        if (c < MAX_CHAR && close_to_open[c] != INVALID_CHAR) {
            const unsigned char open = static_cast<unsigned char>(close_to_open[c]);
            auto& stack = stacks[open];

            if (stack.empty()) {
                THROW_ERROR("Invalid RNA structure: bracket '" + std::string(1, structure[i]) +
                            "' at index " + std::to_string(i) + " was never opened");
            }

            // get's the index of the opening bracket
            const size_t j = stack.back();
            stack.pop_back();

            // record the pairing
            pairings[i] = j;
            pairings[j] = i;

            // Validate pairings if sequence information is available
            if (check_pairs && !can_pair(unmodified_sequence[j], unmodified_sequence[i])) {
                if (have_mod_seq) {
                    std::cerr << "Warning: Base pair '" << mod_sequence[j] << "' can't pair with '"
                              << mod_sequence[i] << "' at indices " << j << ", " << i << '\n';
                } else {
                    std::cerr << "Warning: Base pair '" << unmodified_sequence[j]
                              << "' can't pair with '" << unmodified_sequence[i] << "' at indices "
                              << j << ", " << i << '\n';
                }
            }

            continue;
        }

        // If we reach here, the character is invalid (not a dot or a recognized bracket)
        THROW_ERROR("Invalid RNA structure: invalid character '" + std::string(1, structure[i]) +
                    "' at index " + std::to_string(i));
    }

    // After processing, all stacks should be empty if the structure is well-formed
    for (size_t c = 0; c < MAX_CHAR; ++c) {
        if (!stacks[c].empty()) {
            THROW_ERROR("Invalid RNA structure: opening bracket '" +
                        std::string(1, static_cast<char>(c)) + "' at index " +
                        std::to_string(stacks[c].back()) + " was not closed");
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
// Traverses from right to left to ensure it's sorted by start index without needing an explicit
// sort step.
std::vector<ClosedRegion> RNAProcessor::compute_closed_regions(
    const std::vector<size_t>& pairings) {
    std::vector<ClosedRegion> closed_regions;
    std::stack<ClosedRegion> stack;
    const size_t n = pairings.size();

    for (size_t i = n; i-- > 0;) {
        size_t paired_idx = pairings[i];
        if (paired_idx == NULL_INDEX) continue;

        if (i > paired_idx) {
            // Closing base: push raw pair; its left boundary may be extended later.
            stack.emplace(paired_idx, i);
            continue;
        }

        const ClosedRegion current_pair(i, paired_idx);
        size_t smallest_left = i;

        // Merge any nested regions whose end lies within [i, close)
        while (!stack.empty() && (stack.top().end < current_pair.end)) {
            smallest_left = std::min(smallest_left, stack.top().begin);
            std::cout << "Merging nested region [" << stack.top().begin << ", " << stack.top().end
                      << "] into current pair [" << current_pair.begin << ", " << current_pair.end
                      << "]\n";
            stack.pop();
        }

        if (stack.empty()) THROW_ERROR("Unbalanced pairings\n");

        // Extend region if any nested intervals were merged.
        stack.top().begin = std::min(smallest_left, stack.top().begin);

        // If this opening base starts the region, move it to result.
        if (current_pair.begin == stack.top().begin) {
            closed_regions.push_back(stack.top());
            stack.pop();
        }
    }

    // Because we scanned right-to-left, completed regions were appended in descending begin order.
    // Reverse is linear, not sorting.
    std::reverse(closed_regions.begin(), closed_regions.end());
    // print
    for (const ClosedRegion& cr : closed_regions) {
        std::cerr << "Closed region: [" << cr.begin << ", " << cr.end << "]\n";
    }
    return closed_regions;
}

// AAAAAAAAAAUUUUAAUUUUUU
// (...(...[..)..{..]..)}

// ([...)] = 5, 6, -1, -1, -1, 0, 1
std::vector<size_t> RNAProcessor::compute_cr_pairings(
    const std::vector<ClosedRegion>& closed_regions, size_t rna_size) {
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
    const std::vector<std::string_view>& modified_sequence_views, const all_mod_params& mp,
    const size_t rna_length, bool& has_modified_bases) {
    std::string unmodified_sequence;
    unmodified_sequence.reserve(rna_length);

    bool has_modified_bases_local = false;

    // Convert modified sequence to unmodified sequence
    for (const std::string_view& mod_base : modified_sequence_views) {
        if (mod_base.empty()) {
            THROW_ERROR(
                "A string view in the modified sequence is empty. This should never happen.");
        }

        // If it's an unmodified base, just append it
        if (is_unmod_base(mod_base)) {
            unmodified_sequence.append(mod_base);
            continue;
        }

        has_modified_bases_local = true;

        // If it's a modified base, convert to unmodified using parameters
        const std::string* unmod_base_ptr = mp.get_unmodified_base(std::string(mod_base));

        // lookup modified -> unmodified
        if (unmod_base_ptr != nullptr) {
            unmodified_sequence += *unmod_base_ptr;
        } else {
            THROW_ERROR("Base '" + std::string(mod_base) + "' at index " +
                        std::to_string(unmodified_sequence.size()) +
                        " is not a valid base (not found in modified params and is not a standard "
                        "unmodified base)");
        }
    }

    has_modified_bases = has_modified_bases_local;
    return unmodified_sequence;
}

// Check if a base is unmodified (using lookup table for speed)
bool RNAProcessor::is_unmod_base(const std::string_view& b) {
    if (b.size() != 1)
        return false;  // Multi-character "base" are considered modified (e.g. ❤️‍🩹)
    return unmod_lookup[(unsigned char) b[0]] != 0;
}

bool RNAProcessor::is_unmod_base(char b) {
    return unmod_lookup[static_cast<unsigned char>(b)] != 0;
}

}  // namespace knotergy