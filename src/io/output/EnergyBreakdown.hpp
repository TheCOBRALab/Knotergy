#pragma once

#include "io/output/colors.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

#include <sstream>

namespace knotergy {
class EnergyBreakdown {
   public:
    static std::string node_energy_breakdown(const LoopNode* node,
                                             const ProcessedRNAEntry& rna_entry) {
        size_t max_idx = rna_entry.size();
        LoopType loop_type = node->loop_type;
        bool is_inf = node->is_inf;
        double energy = node->energy;
        size_t begin = node->begin;
        size_t end = node->end;

        std::ostringstream out;
        const bool use_color = should_use_color();

        const auto color = [use_color](const char* ansi_code) -> const char* {
            return use_color ? ansi_code : "";
        };

        const auto idx_w = static_cast<unsigned short>(std::to_string(max_idx).size());

        // Fixed width of: "(<idx_w>, <idx_w>) "
        const auto range_w = static_cast<unsigned short>(1 + idx_w + 2 + idx_w + 2);

        // Loop name
        out << color(ANSI_COLOR_CYAN) << std::left << loop_name(loop_type)
            << color(ANSI_COLOR_RESET) << ' ';

        // Range column
        out << color(ANSI_COLOR_BRIGHT);

        if (loop_type == LoopType::External) {
            out << std::string(range_w, ' ');
        } else {
            out << '(' << std::right << std::setw(idx_w) << begin << ", " << std::right
                << std::setw(idx_w) << end << ") ";
        }

        out << color(ANSI_COLOR_RESET);

        // Energy column
        out << ": " << color(ANSI_COLOR_GREEN);

        if (is_inf) {
            out << std::right << std::setw(9) << "INF";
        } else {
            out << std::right << std::setw(9) << std::fixed << std::setprecision(2) << energy;
        }

        out << color(ANSI_COLOR_RESET) << '\n';

        return out.str();
    }
};

}  // namespace knotergy
