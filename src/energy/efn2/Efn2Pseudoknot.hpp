#pragma once

#include "loop_tree/LoopNode.hpp"

namespace knotergy {

/**
 * Count polymer features directly from a pseudoknot node built by LoopFactory.
 * Use the node's unpaired count, band pairs and children. Stacks and single-base
 * bulges form one helix; other gaps or nested branches split it. Count each band
 * segment once, and skip ordinary nested interiors.
 *
 * Pseudoknotted children need a maximum-cardinality noncrossing choice of their
 * bands to distinguish nested elements from removed pairs. The tree is unchanged.
 */
[[nodiscard]] Efn2PkFeatures get_efn2_features(const LoopNode& node);

/**
 * Hardcoded ShapeKnots/efn2 polymer penalty in 0.01 kcal/mol units.
 * Each closed region is evaluated independently. This intentionally does not
 * reproduce RNAstructure 6.6's accumulated spanning sum across separate knots.
 * Nearest-neighbor band energies are evaluated separately.
 */
[[nodiscard]] double efn2_pk_penalty(const Efn2PkFeatures& features);

}  // namespace knotergy
