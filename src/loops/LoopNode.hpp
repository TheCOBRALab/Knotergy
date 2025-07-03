#pragma once

#include <algorithm>
#include <memory>

#include "../preprocessing/Band.hpp"
#include "../preprocessing/RNAEntry.hpp"

namespace knotergy {
enum class LoopType { Stack, Hairpin, Internal, Multibranch, External, Pseudoknot };
enum class PseudoNestedType { None, InsideBand, OutsideBand, InsideMultiloop };

struct LoopNode {
   public:
    LoopNode(ClosedRegion cr) : begin{cr.begin}, end{cr.end} {}
    LoopNode() : begin{NULL_INDEX}, end{NULL_INDEX} {}

    size_t begin;
    size_t end;

    LoopType loop_type;
    PseudoNestedType pseudo_type = PseudoNestedType::None;
    int number_of_unpaired_bases = 0;
    int number_of_children_inside_band = 0;
    int number_of_children_outside_band = 0;
    int number_of_unpaired_bases_in_children_outside_band = 0;

    std::weak_ptr<LoopNode> parent;
    std::vector<std::shared_ptr<LoopNode>> children;
    std::vector<Band> bands;
    int number_of_bands;
};

inline std::ostream& operator<<(std::ostream& os, const LoopNode& node) {
    os << "LoopNode {\n";
    os << "  begin: " << node.begin << "\n";
    os << "  end: " << node.end << "\n";
    os << "  loop_type: ";

    switch (node.loop_type) {
        case LoopType::Stack:
            os << "Stack";
            break;
        case LoopType::Hairpin:
            os << "Hairpin";
            break;
        case LoopType::Internal:
            os << "Internal";
            break;
        case LoopType::Multibranch:
            os << "Multi";
            break;
        case LoopType::External:
            os << "External";
            break;
        case LoopType::Pseudoknot:
            os << "Pseudoknot";
            break;
    }
    os << "\n";

    os << "  pseudo_type: ";
    switch (node.pseudo_type) {
        case PseudoNestedType::None:
            os << "None";
            break;
        case PseudoNestedType::InsideBand:
            os << "InsideBand";
            break;
        case PseudoNestedType::OutsideBand:
            os << "OutsideBand";
            break;
        case PseudoNestedType::InsideMultiloop:
            os << "InsideMultiloop";
            break;
    }
    os << "\n";

    os << "  number_of_unpaired_bases: " << node.number_of_unpaired_bases << "\n";
    os << "  number_of_children_inside_band: " << node.number_of_children_inside_band << "\n";
    os << "  number_of_children_outside_band: " << node.number_of_children_outside_band << "\n";
    os << "  number_of_unpaired_bases_in_children_outside_band: "
       << node.number_of_unpaired_bases_in_children_outside_band << "\n";

    os << "  number_of_bands: " << node.number_of_bands << "\n";
    os << "  bands:\n";
    for (const auto& band : node.bands) {
        os << "    Band(" << band.left_border << ", " << band.left_inner << ", " << band.right_inner
           << ", " << band.right_border << ")\n";
    }

    os << "  children count: " << node.children.size() << "\n";
    for (const auto& child : node.children) {
        os << "    Child -> begin: " << child->begin << ", end: " << child->end << "\n";
    }

    os << "}\n";
    return os;
}

}  // namespace knotergy