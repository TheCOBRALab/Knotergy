#include "Efn2Pseudoknot.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace knotergy {
namespace {

bool continues_helix(const PKBasePair& bp) {
    if (bp.next == nullptr || !bp.children.empty()) return false;
    const PKBasePair& next = *bp.next;
    const std::size_t left_gap = next.i - bp.i - 1;
    const std::size_t right_gap = bp.j - next.j - 1;
    return left_gap + right_gap < 2;  // Stack or a single-base bulge.
}

using Scores = std::vector<std::size_t>;

bool starts_before(const Band& band, std::size_t position) {
    return band.left_border() < position;
}

std::size_t first_band(const std::vector<Band>& bands, std::size_t position) {
    return std::lower_bound(bands.begin(), bands.end(), position, starts_before)
           - bands.begin();
}

// Maximum weight from 'start' onward, inside the given right boundary.
// 'inside' contains the previously computed optimal nested weights.
std::size_t solve_region(const std::vector<Band>& bands,
                         const Scores& after, const Scores& inside,
                         Scores& dp, std::size_t start, std::size_t right) {
    const std::size_t stop = first_band(bands, right);
    dp[stop] = 0;

    for (std::size_t i = stop; i-- > start;) {
        dp[i] = dp[i + 1];

        if (bands[i].right_border() <= right) {
            dp[i] = std::max(dp[i],
                bands[i].base_pairs().size() + inside[i] + dp[after[i]]);
        }
    }
    return stop;
}

std::vector<const Band*> noncrossing_bands(const LoopNode& node) {
    const auto& bands = node.bands;
    const std::size_t n = bands.size();

    Scores after(n), inside(n), dp(n + 1);

    for (std::size_t i = 0; i < n; ++i) {
        after[i] = first_band(bands, bands[i].right_border());
    }

    // Compute optimal weight nested inside each band.
    for (std::size_t i = n; i-- > 0;) {
        solve_region(bands, after, inside, dp, i + 1, bands[i].right_border());
        inside[i] = dp[i + 1];
    }

    // Recover the optimal noncrossing bands.
    std::vector<const Band*> kept;
    std::vector<std::pair<std::size_t, std::size_t>> pending{{0, node.end}};

    while (!pending.empty()) {
        const auto [start, right] = pending.back();
        pending.pop_back();

        const std::size_t stop = solve_region(bands, after, inside, dp, start, right);

        for (std::size_t i = start; i < stop; ++i) {
            if (bands[i].right_border() <= right &&
                bands[i].base_pairs().size() + inside[i] + dp[after[i]] == dp[i]) {
                kept.push_back(&bands[i]);
                pending.emplace_back(after[i], right);
                pending.emplace_back(i + 1, bands[i].right_border());
                break;
            }
        }
    }
    return kept;
}

using Interval = std::pair<std::size_t, std::size_t>;

void collect_nested_elements(const LoopNode& node, std::vector<Interval>& elements) {
    if (node.loop_type != LoopType::Pseudoknot) {
        elements.emplace_back(node.begin, node.end);
        return;
    }

    for (const Band* band : noncrossing_bands(node))
        elements.emplace_back(band->left_border(), band->right_border());

    for (const LoopNode* child : node.children)
        collect_nested_elements(*child, elements);
}

void count_nested_pseudoknot(const LoopNode& child, Efn2PkFeatures& features) {
    std::vector<Interval> elements;
    collect_nested_elements(child, elements);
    std::sort(elements.begin(), elements.end());

    features.ss += static_cast<int>(child.end - child.begin + 1);

    std::size_t position = child.begin;
    for (const auto& [begin, end] : elements) {
        if (begin < position) continue;

        features.ss -= static_cast<int>(end - begin + 1);
        ++features.ne;
        position = end + 1;
    }
}
}  // namespace

Efn2PkFeatures get_efn2_features(const LoopNode& node) {
    Efn2PkFeatures features;
    features.ss = node.exclusive_unpaired_bases_count;

    for (const Band& band : node.bands) {
        std::size_t length = 1;
        for (const PKBasePair& bp : band.base_pairs()) {
            if (continues_helix(bp)) {
                if (!bp.is_stack(*bp.next)) --features.ss;  // Absorb a single-base bulge.
                ++length;
            } else {
                ++features.il[std::min(length, std::size_t{15})];
                length = 1;
            }
        }
    }

    for (const LoopNode* child : node.children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            count_nested_pseudoknot(*child, features);
        } else {
            ++features.ne;
        }
    }
    return features;
}

double efn2_pk_penalty(const Efn2PkFeatures& features) {
    constexpr double p1 = 0.35;
    constexpr double p2 = 0.65;
    constexpr double a = 6.5;
    constexpr double b = 15.0;
    // RNAstructure 6.6 data_tables/pseudconst.dat; entries 0 and 1 are unused here.
    constexpr std::array<double, 16> c = {0.0,    0.0,    23.9,    23.9,   34.0,  23.8,
                                          34.6,   57.6,   82.1,    65.3,   527.5, 2446.9,
                                          4198.8, 6564.3, 12540.0, 30000.0};

    double span = 0.0;
    for (std::size_t n = 2; n < c.size(); ++n) span += c[n] * c[n] * features.il[n];
    if (span == 0.0) span = 1.0;
    const double contour = features.ss == 0 && features.ne == 0
                               ? p1
                               : p1 * std::log(a * a * features.ss + b * b * features.ne);
    return 100.0 * (contour + p2 * std::log(span));
}

}  // namespace knotergy
