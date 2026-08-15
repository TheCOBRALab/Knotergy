#include "ModInternal.hpp"

#include "energy/modified_bases/ModStack.hpp"

// All the commented stuff is things that should be implemented but isn't currently
// implemented in ViennaRNA. Our goal is to align with ViennaRNA's energy model,
// so we will not implement these until they are implemented in ViennaRNA.

namespace knotergy {

int ModInternal::find_mod_internal_energy(size_t i, size_t j, size_t ci, size_t cj,
                                          const std::vector<std::string_view>& mod_sequence,
                                          const std::string& sequence, vrna_md_param& vp,
                                          const all_mod_params& mp) {
    int unmod_energy = ViennaFunctions::internal_loop_energy(i, j, ci, cj, sequence, vp);
    std::vector<std::string_view> closing_unique_mod_bases =
        ModBaseUtils::unique_modified_bases_at_inner_edge(i, j, mod_sequence);
    std::vector<std::string_view> nested_unique_mod_bases =
        ModBaseUtils::unique_modified_bases_at_inner_edge(ci, cj, mod_sequence);
    if (closing_unique_mod_bases.empty() && nested_unique_mod_bases.empty()) {
        return unmod_energy;
    }

    int n1 = static_cast<int>(ci - i - 1);  // unpaired bases 5' side
    int n2 = static_cast<int>(j - cj - 1);  // unpaired bases 3' side
    int nl = std::max(n1, n2);
    int ns = std::min(n1, n2);
    if (nl == 0) {
        // This is a stack, so we can use the modified stack energy function
        return ModStack::find_mod_stack_energy(i, j, ci, cj, sequence, mod_sequence, vp, mp);
    }

    int e = unmod_energy;
    unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
    unsigned int type2 = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj], vp.md);

    switch (ns) {
        case 0:
            if (nl == 1) {
                // // This is a bulge of size 1, so we can use the modified stack energy function
                // int unmod_stack_energy = ViennaFunctions::stack_energy(i, j, ci, cj, sequence,
                // vp); int mod_stack_energy   = ModStack::find_mod_stack_energy(i, j, ci, cj,
                // sequence, mod_sequence, vp, mp); e += mod_stack_energy - unmod_stack_energy;
            } else {
                if (type1 > 2) {
                    std::string AU_key = ModBaseUtils::join_string_views({i, j}, mod_sequence);
                    int mod_AU_penalty =
                        ModBaseUtils::get_mod_energy(AU_key, closing_unique_mod_bases, mp,
                                                     vp.p->TerminalAU, ModLookup::TerminalAU);
                    e += mod_AU_penalty - vp.p->TerminalAU;
                }
                if (type2 > 2) {
                    std::string AU_key = ModBaseUtils::join_string_views({ci, cj}, mod_sequence);
                    int mod_AU_penalty =
                        ModBaseUtils::get_mod_energy(AU_key, nested_unique_mod_bases, mp,
                                                     vp.p->TerminalAU, ModLookup::TerminalAU);
                    e += mod_AU_penalty - vp.p->TerminalAU;
                }
            }
            break;
            // case 1:
            //     if (nl > 2) {
            //         // Handle the case where the internal loop is a bulge of size 1
            //         e += internal_mismatch_diff(i, j, ci, cj, MismatchType::I1n, mod_sequence,
            //         sequence, vp, mp);
            //     }
            //     break;
            // case 2:
            //     if (nl == 3) {
            //         // Handle the case where the internal loop is a 2x3 internal loop
            //         e += internal_mismatch_diff(i, j, ci, cj,MismatchType::I23,
            //         mod_sequence,sequence, vp, mp);
            //     }
            //     /* fall through */
            // default:
            //     // Handle the case where the internal loop is a general internal loop
            //     e += internal_mismatch_diff(i, j, ci, cj, MismatchType::I, mod_sequence,
            //     sequence, vp, mp); break;
    }
    return e;
}

// enum class MismatchType {I1n, I23, I};
// int internal_mismatch_diff(size_t i, size_t j, size_t ci, size_t cj, MismatchType mm_type,
//                   const std::vector<std::string_view>& mod_sequence,
//                   const std::string& sequence, vrna_md_param& vp,
//                   const all_mod_params& mp) {
//     unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
//     unsigned int type2 = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj], vp.md);
//     int si1 = ViennaUtils::fast_nucleotide_encode(sequence[i + 1]);   // 5' mismatch nt
//     of closing pair int sj1 = ViennaUtils::fast_nucleotide_encode(sequence[j - 1]);   //
//     3' mismatch nt of closing pair int sp1 = ViennaUtils::fast_nucleotide_encode(sequence[ci -
//     1]);  // 5' mismatch nt of enclosed pair int sq1 =
//     ViennaUtils::fast_nucleotide_encode(sequence[cj + 1]);  // 3' mismatch nt of enclosed
//     pair

//     int unmod_mismatch1;
//     int unmod_mismatch2;

//     switch (mm_type) {
//         case MismatchType::I1n:
//             unmod_mismatch1 = vp.p->mismatch1nI[type1][si1][sj1];
//             unmod_mismatch2 = vp.p->mismatch1nI[type2][sq1][sp1];
//             break;

//         case MismatchType::I23:
//             unmod_mismatch1 = vp.p->mismatch23I[type1][si1][sj1];
//             unmod_mismatch2 = vp.p->mismatch23I[type2][sq1][sp1];
//             break;

//         case MismatchType::I:
//             unmod_mismatch1 = vp.p->mismatchI[type1][si1][sj1];
//             unmod_mismatch2 = vp.p->mismatchI[type2][sq1][sp1];
//             break;

//     }

//     std::string mismatch1_key = ModBaseUtils::join_string_views({i, j, i + 1}, mod_sequence);
//     std::string mismatch2_key = ModBaseUtils::join_string_views({ci, cj, ci - 1}, mod_sequence);
//     std::vector<std::string_view> closing_unique_mod_bases =
//     ModBaseUtils::unique_modified_bases_at_inner_edge(i, j, mod_sequence);
//     std::vector<std::string_view> nested_unique_mod_bases  =
//     ModBaseUtils::unique_modified_bases_at_inner_edge(ci, cj, mod_sequence); int mod_mismatch1 =
//     ModBaseUtils::get_mod_energy(mismatch1_key, closing_unique_mod_bases, mp, unmod_mismatch1,
//     ModLookup::Mismatch); int mod_mismatch2 = ModBaseUtils::get_mod_energy(mismatch2_key,
//     nested_unique_mod_bases, mp, unmod_mismatch2, ModLookup::Mismatch);

//     return (mod_mismatch1 + mod_mismatch2) - (unmod_mismatch1 + unmod_mismatch2);
// }

int ModInternal::find_mod_internal_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                          const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                          const all_mod_params& mp) {
    return find_mod_internal_energy(bp.i, bp.j, next_bp.i, next_bp.j, pRNA.get_modified_sequence(),
                                    pRNA.get_sequence(), vp, mp);
}

}  // namespace knotergy