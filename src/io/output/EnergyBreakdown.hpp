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

enum class VerbosityLevel { Quiet = 0, Verbose = 1, Detailed = 2 };

class EnergyBreakdown {
   public:
    static std::string node_energy_breakdown(const LoopNode* node,
                                             const ProcessedRNAEntry& rna_entry,
                                             VerbosityLevel verbosity = VerbosityLevel::Quiet) {
        if (verbosity == VerbosityLevel::Quiet) {
            return "";
        }

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

        // ViennaRNA aligns ':' at column 45 (44 characters before it).
        constexpr std::size_t DESCRIPTION_WIDTH = 44;

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

        // ------------------------------------------------------------
        // Detailed pseudoknot band breakdown
        // ------------------------------------------------------------
        if (node->loop_type == LoopType::Pseudoknot && verbosity == VerbosityLevel::Detailed) {
            const std::string pk_label = "Pseudoknot-level energy (excludes bands)";

            // Longest band row:
            //
            //     Interior loop (  86, 886) AU; (  88, 884) GU
            //
            // 4  = indentation
            // 13 = loop_name() width
            // first pair  = 7 + 2 * idx_w
            // second pair = 8 + 2 * idx_w
            const std::size_t band_description_width = 4 + 13 + (7 + 2 * idx_w) + (8 + 2 * idx_w);

            // --------------------------------------------------------
            // Pseudoknot-level energy
            // --------------------------------------------------------
            out << "    " << color(ANSI_COLOR_BLUE) << pk_label << color(ANSI_COLOR_RESET);

            const std::size_t pk_visible_width = 4 + pk_label.size();

            if (pk_visible_width < band_description_width) {
                out << std::string(band_description_width - pk_visible_width, ' ');
            }

            out << ": " << color(ANSI_COLOR_GRAY) << std::right << std::setw(5)
                << std::llround(node->pk_level_energy) << color(ANSI_COLOR_RESET) << '\n';

            // --------------------------------------------------------
            // Bands
            // --------------------------------------------------------
            for (std::size_t band_idx = 0; band_idx < node->bands.size(); ++band_idx) {
                const Band& band = node->bands[band_idx];

                out << "  " << color(ANSI_COLOR_BLUE) << "Band " << band_idx + 1
                    << color(ANSI_COLOR_RESET) << " (" << band.left_border() + 1 << ", "
                    << band.left_inner() + 1 << ", " << band.right_inner() + 1 << ", "
                    << band.right_border() + 1 << ")\n";

                const std::vector<PKBasePair>& bps = band.base_pairs();

                // ----------------------------------------------------
                // Base pairs in this band
                // ----------------------------------------------------
                for (std::size_t idx = 0; idx < bps.size(); ++idx) {
                    const PKBasePair& bp = bps[idx];

                    std::ostringstream bp_details;
                    std::ostringstream styled_bp_details;

                    // Current base pair
                    bp_details << " (" << std::right << std::setw(static_cast<int>(idx_w))
                               << bp.i + 1 << "," << std::setw(static_cast<int>(idx_w)) << bp.j + 1
                               << ") " << sequence[bp.i] << sequence[bp.j];

                    styled_bp_details << " (" << std::right << std::setw(static_cast<int>(idx_w))
                                      << bp.i + 1 << "," << std::setw(static_cast<int>(idx_w))
                                      << bp.j + 1 << ") " << color(ANSI_COLOR_BRIGHT)
                                      << sequence[bp.i] << sequence[bp.j]
                                      << color(ANSI_COLOR_RESET);

                    // Stack and internal loops are defined relative to
                    // the next base pair in the band.
                    if ((bp.loop_type == LoopType::Stack || bp.loop_type == LoopType::Internal) &&
                        idx + 1 < bps.size()) {
                        const PKBasePair& next_bp = bps[idx + 1];

                        bp_details << "; (" << std::setw(static_cast<int>(idx_w)) << next_bp.i + 1
                                   << "," << std::setw(static_cast<int>(idx_w)) << next_bp.j + 1
                                   << ") " << sequence[next_bp.i] << sequence[next_bp.j];

                        styled_bp_details
                            << "; (" << std::setw(static_cast<int>(idx_w)) << next_bp.i + 1 << ","
                            << std::setw(static_cast<int>(idx_w)) << next_bp.j + 1 << ") "
                            << color(ANSI_COLOR_BRIGHT) << sequence[next_bp.i]
                            << sequence[next_bp.j] << color(ANSI_COLOR_RESET);
                    }

                    const std::size_t bp_visible_width =
                        4 + std::string(loop_name(bp.loop_type)).size() + bp_details.str().size();

                    out << "    " << color(ANSI_COLOR_BLUE) << loop_name(bp.loop_type)
                        << color(ANSI_COLOR_RESET) << styled_bp_details.str();

                    if (bp_visible_width < band_description_width) {
                        out << std::string(band_description_width - bp_visible_width, ' ');
                    }

                    out << ": " << color(ANSI_COLOR_GRAY) << std::right << std::setw(5)
                        << std::llround(bp.energy) << color(ANSI_COLOR_RESET) << '\n';
                }
            }
        }

        return out.str();
    }
};

}  // namespace knotergy