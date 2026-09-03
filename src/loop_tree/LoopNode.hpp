#pragma once

#include "energy/pseudoknots/PKEnergyBreakdown.hpp"
#include "io/output/colors.hpp"
#include "loop_tree/LoopTypes.hpp"
#include "loop_tree/bands/Band.hpp"
#include "preprocessing/ClosedRegion.hpp"

#include <algorithm>
#include <iomanip>
#include <memory>
#include <sstream>

namespace knotergy {

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
    //  -------------- Constructors -------------------
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

    //  ------------ Member Variables -------------------
    size_t begin;  ///< 5' boundary position (or NULL_INDEX for external loop).
    size_t end;    ///< 3' boundary position (or NULL_INDEX for external loop).

    LoopType loop_type = LoopType::Unknown;                 ///< Type of this loop.
    PseudoNestedType pseudo_type = PseudoNestedType::None;  ///< Pseudoknot nesting type.
    int exclusive_unpaired_bases_count = 0;                 ///< Unpaired bases only in this loop.
    int total_unpaired_bases_count = 0;      ///< Unpaired bases in loop + nested children.
    int number_of_withinband_children = 0;   ///< Count of children within pseudoknot bands.
    int number_of_outsideband_children = 0;  ///< Count of nested children.
    LoopNode* parent = nullptr;              ///< Parent loop node (weak to avoid cycles).
    std::vector<LoopNode*> children;         ///< Child loop nodes. (sorted by start)
    std::vector<Band> bands;                 ///< Pseudoknot bands (empty if not pseudoknotted).
    int total_number_of_base_pairs = 1;  ///< # of pairs in this closed region (excluding children)
    double energy = 0;                   ///< Computed energy (set by ComputeEnergy).
    bool is_inf = false;                 ///< Flag for infinite energy (e.g., invalid structures).
    PKEnergyBreakdown pk_energy_breakdown;  ///< Detailed breakdown of pseudoknot energy.
};

/**
 * @brief Stream insertion operator for LoopNode.
 *
 * Prints detailed information about the loop node including type, boundaries,
 * unpaired base counts, bands, and children.
 *
 * How to use:
 * LoopNode node = ...;
 * std::cout << node;
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
    os << pk_nested_type(node.pseudo_type);
    os << "\n";

    os << "  exclusive_unpaired_bases_count: " << node.exclusive_unpaired_bases_count << "\n";
    os << "  total_unpaired_bases_count: " << node.total_unpaired_bases_count << "\n";
    os << "  number_of_children_inside_band: " << node.number_of_withinband_children << "\n";
    os << "  number_of_outsideband_children: " << node.number_of_outsideband_children << "\n";
    os << "  number_of_bands: " << node.bands.size() << "\n";
    os << "  bands:\n";
    for (const Band& band : node.bands) {
        os << "    Band(" << band.left_border() << ", " << band.left_inner() << ", "
           << band.right_inner() << ", " << band.right_border() << ")\n";
        for (const PKBasePair& base_pair : band.base_pairs()) {
            os << "        PKBasePair(" << base_pair.i << ", " << base_pair.j << ")\n";
        }
    }
    os << "  children count: " << node.children.size() << "\n";
    for (const LoopNode* child : node.children) {
        os << "    Child -> begin: " << child->begin << ", end: " << child->end << "\n";
    }
    os << "}\n";
    return os;
}

}  // namespace knotergy