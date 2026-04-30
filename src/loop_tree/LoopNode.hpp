#pragma once

#include <algorithm>
#include <iomanip>
#include <memory>
#include <sstream>

#include "../loop_tree/Band.hpp"
#include "../preprocessing/ClosedRegion.hpp"

namespace knotergy {
/**
 * @brief Enumeration of loop types in RNA secondary structures.
 */
enum class LoopType { Unknown, Stack, Hairpin, Internal, Multibranch, External, Pseudoknot };

/**
 * @brief Enumeration of pseudoknot nesting types.
 */
enum class PseudoNestedType { None, WithinBand, Nested };
// Within Band (((..(...).(.[[[.)..)))]]]
// This hairpin     ^   ^ is within a band
// Nested      (((..(...).[[[...)))]]]
// This hairpin     ^   ^ is nested inside a band

/**
 * @brief Get a human-readable name for a loop type.
 *
 * @param t The loop type.
 * @return String representation of the loop type.
 */
static inline const char* loop_name(LoopType t) {
    switch (t) {
        case LoopType::Stack:
            return "Stack    loop";
        case LoopType::Hairpin:
            return "Hairpin  loop";
        case LoopType::Internal:
            return "Internal loop";
        case LoopType::Multibranch:
            return "Multi    loop";
        case LoopType::External:
            return "External loop";
        case LoopType::Pseudoknot:
            return "Pseudo   loop";
        case LoopType::Unknown:
            return "Unknown  loop";
    }
    return "Unknown  loop";
}

/**
 * @brief Represents a node in the loop tree of an RNA secondary structure.
 *
 * Each LoopNode corresponds to a loop region in the RNA structure and contains
 * information about its boundaries, type, nested children, pseudoknot bands,
 * unpaired base counts, and computed energy. The loop tree is built by LoopFactory
 * and processed by ComputeEnergy.
 */
struct LoopNode {
   public:
    /**
     * @brief Construct a LoopNode from a closed region.
     *
     * @param cr Closed region defining the loop boundaries.
     */
    LoopNode(ClosedRegion cr) : begin{cr.begin}, end{cr.end} {}

    /**
     * @brief Default constructor for external loop (no boundaries).
     */
    LoopNode() : begin{NULL_INDEX}, end{NULL_INDEX} {}

    size_t begin;  ///< 5' boundary position (or NULL_INDEX for external loop).
    size_t end;    ///< 3' boundary position (or NULL_INDEX for external loop).

    LoopType loop_type = LoopType::Unknown;                 ///< Type of this loop.
    PseudoNestedType pseudo_type = PseudoNestedType::None;  ///< Pseudoknot nesting type.
    int exclusive_unpaired_bases_count = 0;                 ///< Unpaired bases only in this loop.
    int total_unpaired_bases_count = 0;     ///< Unpaired bases in loop + nested children.
    int number_of_withinband_children = 0;  ///< Count of children within pseudoknot bands.
    int number_of_nested_children = 0;      ///< Count of nested children.
    std::weak_ptr<LoopNode> parent;         ///< Parent loop node (weak to avoid cycles).
    std::vector<std::shared_ptr<LoopNode>> children;  ///< Child loop nodes.
    std::vector<Band> bands;  ///< Pseudoknot bands (empty if not pseudoknotted).
    int number_of_bands = 0;  ///< Number of bands in this loop.
    double energy = 0;        ///< Computed energy (set by ComputeEnergy).
    bool is_inf = false;      ///< Flag for infinite energy (e.g., invalid structures).

    /**
     * @brief Generate a formatted string for energy breakdown output.
     *
     * Used for verbose output showing loop type, position range, and energy.
     *
     * @param max_idx Maximum index in the structure (for formatting width).
     * @return Formatted string with loop information and energy.
     */
    std::string energy_breakdown(size_t max_idx) const {
        std::ostringstream out;

        const unsigned short idx_w = static_cast<unsigned short>(std::to_string(max_idx).size());

        // This is the *fixed* width of: "[<idx_w>, <idx_w>] "
        const unsigned short range_w =
            static_cast<unsigned short>(1 + idx_w + 2 + idx_w + 2);  // '[' + a + ", " + b + "] "

        // Colored name padded
        out << "\x1b[36m" << std::left << loop_name(loop_type) << "\x1b[0m ";

        // Range column (i, j)
        if (loop_type == LoopType::External) {
            out << std::string(range_w, ' ');
        } else {
            out << "(" << std::right << std::setw(idx_w) << begin << ", " << std::right
                << std::setw(idx_w) << end << ") ";
        }

        // Energy column
        if (is_inf) {
            out << ": " << std::right << std::setw(9) << "inf\n";
        } else {
            out << ": " << std::right << std::setw(9) << std::fixed << std::setprecision(2)
                << energy << "\n";
        }

        return out.str();
    }
};

/**
 * @brief Stream insertion operator for LoopNode.
 *
 * Prints detailed information about the loop node including type, boundaries,
 * unpaired base counts, bands, and children.
 *
 * @param os Output stream.
 * @param node LoopNode to print.
 * @return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, const LoopNode& node) {
    os << "LoopNode {\n";
    os << "  begin: " << node.begin << "\n";
    os << "  end: " << node.end << "\n";
    os << "  loop_type: ";
    os << loop_name(node.loop_type);
    os << "\n";

    os << "  pseudo_type: ";
    switch (node.pseudo_type) {
        case PseudoNestedType::None:
            os << "None";
            break;
        case PseudoNestedType::WithinBand:
            os << "WithinBand";
            break;
        case PseudoNestedType::Nested:
            os << "Nested";
            break;
    }
    os << "\n";

    os << "  exclusive_unpaired_bases_count: " << node.exclusive_unpaired_bases_count << "\n";
    os << "  total_unpaired_bases_count: " << node.total_unpaired_bases_count << "\n";
    os << "  number_of_children_inside_band: " << node.number_of_withinband_children << "\n";
    os << "  number_of_nested_children: " << node.number_of_nested_children << "\n";
    os << "  number_of_bands: " << node.number_of_bands << "\n";
    os << "  bands:\n";
    for (const auto& band : node.bands) {
        os << "    Band(" << band.left_border() << ", " << band.left_inner() << ", "
           << band.right_inner() << ", " << band.right_border() << ")\n";
        for (const auto& base_pair : band.base_pairs()) {
            os << "        BasePair(" << base_pair.i << ", " << base_pair.j << ")\n";
        }
    }
    os << "  children count: " << node.children.size() << "\n";
    for (const auto& child : node.children) {
        os << "    Child -> begin: " << child->begin << ", end: " << child->end << "\n";
    }
    os << "}\n";
    return os;
}

}  // namespace knotergy