#include "ViennaDangles.hpp"

namespace knotergy {

enum TouchingRight {
    RightFree = 0,
    RightTaken = 1
};

int ViennaDangles::get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence, vrna_md_t& md) {
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, sequence, md);
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    return process_chains(dangle_chains, children, dangle_energies);
}

int ViennaDangles::get_multi_dangle_1(const LoopNode& node, const std::string& sequence, vrna_md_t& md){
    const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
    bool is_external = false;
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, sequence, md, is_external);
    DangleSet ml_dangle_energy = get_ml_dangle_energy(node, sequence, md);
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    return process_ml_chains(dangle_chains, children, dangle_energies, node, ml_dangle_energy);
}

DangleSet ViennaDangles::get_ml_dangle_energy(const LoopNode& node, const std::string& sequence, vrna_md_t& md){
    vrna_param_t* P = vrna_params(&md);
    size_t pi = node.begin, pj = node.end;
    DangleSet ml_dangle; // closing pair dangles
    int n5d = vrna_nucleotide_encode(sequence[pi + 1], &md);
    int n3d = vrna_nucleotide_encode(sequence[pj - 1], &md);
    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[pi], sequence[pj], md);
    unsigned int rev_pair_type = ViennaUtils::reverse_pair_type(pair_type, md);

    ml_dangle.no_dangle    = vrna_E_multibranch_stem(rev_pair_type, -1,  -1,  P);
    ml_dangle.left_dangle  = vrna_E_multibranch_stem(rev_pair_type, n3d, -1,  P);
    ml_dangle.right_dangle = vrna_E_multibranch_stem(rev_pair_type, -1,  n5d, P);
    ml_dangle.both_dangle  = vrna_E_multibranch_stem(rev_pair_type, n3d, n5d, P);
    if (P) free(P);
    return ml_dangle;
}

std::vector<DangleSet> ViennaDangles::populate_children_dangle_energies(
    const std::vector<std::shared_ptr<LoopNode>>& children,
    const std::string& sequence,
    vrna_md_t& md,
    const bool& is_external
) {
    vrna_param_t* P = vrna_params(&md);
    std::vector<DangleSet> dangle_energies;
    int (*vrna_E_stem)(unsigned int, int, int, vrna_param_t*); 
    vrna_E_stem = is_external ? vrna_E_exterior_stem : vrna_E_multibranch_stem;
    
    for (size_t i = 0; i < children.size(); ++i) {
        const std::shared_ptr<LoopNode>& child = children[i];
        const size_t& ci = child->begin, cj = child->end;

        // Compute all dangle energies for this child
        int n5d = ci > 0 ? vrna_nucleotide_encode(sequence[ci - 1], &md) : -1;
        int n3d = cj < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[cj + 1], &md) : -1;
        unsigned int pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], md);

        DangleSet d_energy;
        d_energy.no_dangle  = vrna_E_stem(pair_type, -1,  -1,  P);
        d_energy.left_dangle  = vrna_E_stem(pair_type, n5d, -1,  P);
        d_energy.right_dangle = vrna_E_stem(pair_type, -1,  n3d, P);
        d_energy.both_dangle  = vrna_E_stem(pair_type, n5d, n3d, P);
        dangle_energies.push_back(d_energy);
    }
    if (P) free(P);
    return dangle_energies;
}

std::vector<std::vector<size_t>> ViennaDangles::get_dangle_chains(const std::vector<std::shared_ptr<LoopNode>>& children) {
    std::vector<std::vector<size_t>> dangle_chains;
    dangle_chains.reserve(children.size());
    dangle_chains.push_back({0});
    for (size_t i = 1; i < children.size(); ++i) {
        const std::shared_ptr<LoopNode>& child = children[i];
        const std::shared_ptr<LoopNode>& prev_child = children[i - 1];
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }

        // if there's 1 or no unpaired bases between children, they share or have no dangles
        if (child->begin - prev_child->end <= 2) {
            dangle_chains.back().push_back(i);
        } else{
            dangle_chains.push_back({i});
        }
    }
    return dangle_chains;
}

int ViennaDangles::process_chain(
                    const std::vector<size_t>& chain,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const bool& disable_first_left_dangle,
                    const bool& disable_last_right_dangle,
                    std::array<int,2> prev_init
                ) {
    if (chain.size() == 1) {
        if (disable_first_left_dangle && disable_last_right_dangle) {
            return dangle_energies[chain.front()].no_dangle;
        } else if (disable_first_left_dangle) {
            return std::min(dangle_energies[chain.front()].no_dangle,
                            dangle_energies[chain.front()].right_dangle);
        } else if (disable_last_right_dangle) {
            return std::min(dangle_energies[chain.front()].no_dangle,
                            dangle_energies[chain.front()].left_dangle);
        }
        return dangle_energies[chain.front()].min();
    }
    int dangle_energy = 0;
    std::array<int,2> prev = prev_init; // default {0, INF}
    for (size_t idx : chain) {
        const DangleSet& energies = dangle_energies[idx];
        std::array<int,2> cur = {INF, INF};

        // Check if left or right dangle is possible based on adjacency (no unpaired bases in between)
        bool disable_left_dangle = idx != chain.front() && contiguous_children(*children[idx], *children[idx - 1]) || (idx == chain.front() && disable_first_left_dangle);
        bool disable_right_dangle = idx != chain.back() && contiguous_children(*children[idx], *children[idx + 1]) || (idx == chain.back() && disable_last_right_dangle);

        // energies: no dangle, left dangle, right dangle, both dangles
        int eNone = energies.no_dangle, eLeft = energies.left_dangle, eRight = energies.right_dangle, eBoth = energies.both_dangle;

        // Update touching_right based on previous state and current possibilities
        if (disable_left_dangle && disable_right_dangle) {
            cur[RightFree] = prev[RightFree] + eNone;
            // cur[RightTaken] = INF; // impossible
        } else if (disable_left_dangle) {
            cur[RightFree] = prev[RightFree] + eNone; // impossible case: prev[RightTaken] + eNone
            cur[RightTaken] = std::min({prev[RightFree] + eRight, prev[RightTaken] + eRight});
        } else if (disable_right_dangle) {
            cur[RightFree] = std::min({prev[RightFree] + eNone, prev[RightFree] + eLeft, prev[RightTaken] + eNone, prev[RightTaken] + eLeft});
            // cur[RightTaken] = INF; // impossible
        } else {
            cur[RightFree] = std::min({prev[RightFree] + eNone , prev[RightFree] + eLeft,  prev[RightTaken] + eNone});
            cur[RightTaken] = std::min({prev[RightFree] + eLeft , prev[RightFree] + eRight, prev[RightTaken] + eRight, prev[RightFree] + eBoth});
        }
        
        // Sanity check for overflow
        if ((cur[RightFree] > INF) || (cur[RightTaken] > INF)) {
            THROW_ERROR("Dangle energy overflow detected in external loop dangle calculation.");
        }
        prev = cur;
    }
    return std::min(prev[RightFree], prev[RightTaken]);
}

int ViennaDangles::process_chains(
                    const std::vector<std::vector<size_t>>& dangle_chains,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const bool& disable_first_left_dangle,
                    const bool& disable_last_right_dangle,
                    std::array<int,2> init
                ) {
    int dangle_energy = 0;
    bool first = true;
    for (const std::vector<size_t>& chain : dangle_chains) {
        if (first) {
            first = false;
            dangle_energy += process_chain(chain, children, dangle_energies, disable_first_left_dangle, false, init);
            continue;
        }
        dangle_energy += process_chain(chain, children, dangle_energies);
    }
    return dangle_energy;
}

int ViennaDangles::process_ml_chains(
                    const std::vector<std::vector<size_t>>& dangle_chains,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const LoopNode& node,
                    const DangleSet ml_dangle_energy
                ) {
    enum class ML_Type { None, Front, Back, Both, Loop};
   
    bool front_dangle = children.front()->begin - node.begin <= 2;
    bool back_dangle  = node.end - children.back()->end <= 2;
    int best_left  = std::min(ml_dangle_energy.no_dangle, ml_dangle_energy.left_dangle);
    int best_right = std::min(ml_dangle_energy.no_dangle, ml_dangle_energy.right_dangle);

    ML_Type ml_type;
    if (front_dangle && back_dangle && children.size() == 1){
        ml_type = ML_Type::Loop;
    } else if (front_dangle && back_dangle) {
        ml_type = ML_Type::Both;
    } else if (front_dangle) {
        ml_type = ML_Type::Front;
    } else if (back_dangle) {
        ml_type = ML_Type::Back;
    } else {
        ml_type = ML_Type::None;
    }

    if (ml_type == ML_Type::None){
        return ml_dangle_energy.min() + process_chains(dangle_chains, children, dangle_energies);
    }

    if (ml_type == ML_Type::Front){
        if (contiguous(node.begin, children.front()->begin)){
            return process_chains(dangle_chains, children, dangle_energies, true) + best_left;
        } else {
            return process_chains(dangle_chains, children, dangle_energies, false, false, {ml_dangle_energy.no_dangle, best_left});
        }
    }
    
    if (ml_type == ML_Type::Back){
        if (contiguous(node.end, children.back()->end)){
            return process_chains(dangle_chains, children, dangle_energies, false,  true) + best_right;
        } else {
            return process_chains(dangle_chains, children, dangle_energies, false, false, {ml_dangle_energy.no_dangle, best_right});
        }
    }
    
    int dangle_energy = 0;
    return dangle_energy;
}

}