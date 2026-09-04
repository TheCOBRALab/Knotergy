#pragma once

#include "io/output/colors.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace knotergy {

class EnergyBreakdown {
   public:
    static std::string node_energy_breakdown(const LoopNode* node,
                                             const ProcessedRNAEntry& rna_entry) {
        const bool use_color = should_use_color();

        const auto color = [use_color](const char* ansi_code) -> const char* {
            return use_color ? ansi_code : "";
        };

        const std::string& sequence = rna_entry.get_sequence();

        // ViennaRNA uses at least 3 characters for indices.
        const std::size_t idx_w = std::max<std::size_t>(3, std::to_string(rna_entry.size()).size());

        // ------------------------------------------------------------
        // Build loop details
        // ------------------------------------------------------------
        std::ostringstream details;
        std::ostringstream styled_details;
        std::ostringstream out;

        if (node->loop_type != LoopType::External) {
            // Outer pair
            details << " (" << std::right << std::setw(static_cast<int>(idx_w)) << node->begin + 1
                    << "," << std::setw(static_cast<int>(idx_w)) << node->end + 1 << ") "
                    << sequence[node->begin] << sequence[node->end];

            styled_details << " (" << std::right << std::setw(static_cast<int>(idx_w))
                           << node->begin + 1 << "," << std::setw(static_cast<int>(idx_w))
                           << node->end + 1 << ") " << color(ANSI_COLOR_BRIGHT)
                           << sequence[node->begin] << sequence[node->end]
                           << color(ANSI_COLOR_RESET);

            if ((node->loop_type == LoopType::Stack || node->loop_type == LoopType::Internal) &&
                !node->children.empty()) {
                const LoopNode* child = node->children.front();

                details << "; (" << std::setw(static_cast<int>(idx_w)) << child->begin + 1 << ","
                        << std::setw(static_cast<int>(idx_w)) << child->end + 1 << ") "
                        << sequence[child->begin] << sequence[child->end];

                styled_details << "; (" << std::setw(static_cast<int>(idx_w)) << child->begin + 1
                               << "," << std::setw(static_cast<int>(idx_w)) << child->end + 1
                               << ") " << color(ANSI_COLOR_BRIGHT) << sequence[child->begin]
                               << sequence[child->end] << color(ANSI_COLOR_RESET);
            }
        }

        const std::string detail_string = details.str();

        // loop_name() is always 13 visible characters.
        const std::size_t visible_width =
            std::string(loop_name(node->loop_type)).size() + detail_string.size();

        // ViennaRNA aligns ':' at column 41 (40 characters before it).
        constexpr std::size_t DESCRIPTION_WIDTH = 40;

        // ------------------------------------------------------------
        // Loop description
        // ------------------------------------------------------------
        out << color(ANSI_COLOR_CYAN) << loop_name(node->loop_type) << color(ANSI_COLOR_RESET);

        out << styled_details.str();

        if (visible_width < DESCRIPTION_WIDTH) {
            out << std::string(DESCRIPTION_WIDTH - visible_width, ' ');
        }
        // ------------------------------------------------------------
        // Energy
        // ------------------------------------------------------------
        out << ": " << color(ANSI_COLOR_GREEN);

        if (node->is_inf) {
            out << std::right << std::setw(5) << "INF";
        } else {
            out << std::right << std::setw(5) << std::llround(node->energy);
        }

        out << color(ANSI_COLOR_RESET) << '\n';

        return out.str();
    }
};

}  // namespace knotergy