#pragma once

#include <algorithm>
#include <memory>
#include <sstream>
#include <iomanip>

#include "../loop_tree/Band.hpp"
#include "../preprocessing/ClosedRegion.hpp"

namespace knotergy {
enum class LoopType { Stack, Hairpin, Internal, Multibranch, External, Pseudoknot };
enum class PseudoNestedType { None, WithinBand, Nested};
// Within Band (((..(...).(.[[[.)..)))]]]
// This hairpin     ^   ^ is within a band
// Nested      (((..(...).[[[...)))]]]
// This hairpin     ^   ^ is nested inside a band

static inline const char* loop_name(LoopType t) {
    switch (t) {
        case LoopType::Stack:       return "Stack    loop";
        case LoopType::Hairpin:     return "Hairpin  loop";
        case LoopType::Internal:    return "Internal loop";
        case LoopType::Multibranch: return "Multi    loop";
        case LoopType::External:    return "External loop";
        case LoopType::Pseudoknot:  return "Pseudo   loop";
    }
    return "Unknown  loop";
}

struct LoopNode {
   public:
    LoopNode(ClosedRegion cr) : begin{cr.begin}, end{cr.end} {}
    LoopNode() : begin{NULL_INDEX}, end{NULL_INDEX} {}

    size_t begin;
    size_t end;

    LoopType loop_type;
    PseudoNestedType pseudo_type = PseudoNestedType::None;
    int exclusive_unpaired_bases_count = 0; // unpaired bases in loop only
    int total_unpaired_bases_count = 0;  // unpaired bases in loop + unpaired bases in nested children
    int number_of_withinband_children = 0;
    int number_of_nested_children = 0;
    int number_of_unpaired_bases_in_nested_children = 0;

    std::weak_ptr<LoopNode> parent;
    std::vector<std::shared_ptr<LoopNode>> children;
    std::vector<Band> bands;
    int number_of_bands;
    double energy = 0; // calculated in ComputeEnergy.cpp


    std::string energy_breakdown(size_t max_idx) const {
        std::ostringstream out;

        const unsigned short idx_w  = static_cast<unsigned short>(std::to_string(max_idx).size());

        // This is the *fixed* width of: "[<idx_w>, <idx_w>] "
        const unsigned short range_w = 1 + idx_w + 2 + idx_w + 2; // '[' + a + ", " + b + "] "

        // Colored name padded (setw applies only to loop_name, not escape codes)
        out << "\x1b[36m"
            << std::left << loop_name(loop_type)
            << "\x1b[0m ";

        // Range column
        if (loop_type == LoopType::External) {
            out << std::string(range_w, ' ');
        } else {
            out << "("
                << std::right << std::setw(idx_w) << begin
                << ", "
                << std::right << std::setw(idx_w) << end
                << ") ";
        }

        // Energy column
        out << ": "
            << std::right << std::setw(9) << std::fixed << std::setprecision(2) << energy
            << "\n";

        return out.str();
    }
};



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
    os << "  number_of_unpaired_bases_in_nested_children: "
       << node.number_of_unpaired_bases_in_nested_children << "\n";
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