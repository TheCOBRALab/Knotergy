#include "ViennaDangles.hpp"

namespace knotergy {
enum TouchingRight {N = 0, Y = 1};

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
    return 0;
}

DangleSet ViennaDangles::get_ml_dangle_energy(const LoopNode& node, const std::string& sequence, vrna_md_t& md){
    vrna_param_t* P = vrna_params(&md);
    size_t pi = node.begin, pj = node.end;
    DangleSet ml_dangle; // closing pair dangles
    int n5d = vrna_nucleotide_encode(sequence[pi + 1], &md);
    int n3d = vrna_nucleotide_encode(sequence[pj - 1], &md);
    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[pi], sequence[pj], md);
    unsigned int rev_pair_type = ViennaUtils::reverse_pair_type(pair_type, md);

    ml_dangle.no_dangle  = vrna_E_multibranch_stem(rev_pair_type, -1,  -1,  P);
    ml_dangle.left_dangle  = vrna_E_multibranch_stem(rev_pair_type, -1,  n5d, P);
    ml_dangle.right_dangle = vrna_E_multibranch_stem(rev_pair_type, n3d, -1,  P);
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
                    const bool& disable_first_right_dangle
                ) {
    int dangle_energy = 0;
    std::array<int,2> prev_TR = {0, 0};
    for (size_t idx : chain) {
        const DangleSet& energies = dangle_energies[idx];
        std::array<int,2> touching_right = {INF, INF};

        // Check if left or right dangle is possible based on adjacency (no unpaired bases in between)
        bool no_left_dangle = idx != chain.front() && (children[idx]->begin - children[idx - 1]->end == 1);
        bool no_right_dangle = idx != chain.back() && (children[idx + 1]->begin - children[idx]->end == 1);

        // energies: no dangle, left dangle, right dangle, both dangles
        int eNone = energies.no_dangle, eLeft = energies.left_dangle, eRight = energies.right_dangle, eBoth = energies.both_dangle;

        // Update touching_right based on previous state and current possibilities
        if (no_left_dangle && no_right_dangle) {
            touching_right[N] = prev_TR[N] + eNone;
            // touching_right[1] = INF; // impossible
        } else if (no_left_dangle) {
            touching_right[N] = prev_TR[N] + eNone; // impossible case: prev[1] + eNone
            touching_right[Y] = std::min({prev_TR[N] + eRight, prev_TR[Y] + eRight});
        } else if (no_right_dangle) {
            touching_right[N] = std::min({prev_TR[N] + eNone, prev_TR[N] + eLeft, prev_TR[Y] + eNone, prev_TR[Y] + eLeft});
            // touching_right[1] = INF; // impossible
        } else {
            touching_right[N] = std::min({prev_TR[N] + eNone , prev_TR[N] + eLeft,  prev_TR[Y] + eNone});
            touching_right[Y] = std::min({prev_TR[N] + eLeft , prev_TR[N] + eRight, prev_TR[Y] + eRight, prev_TR[N] + eBoth});
        }
        
        // Sanity check for overflow
        if ((touching_right[N] > INF) or (touching_right[Y] > INF)){
            THROW_ERROR("Dangle energy overflow detected in external loop dangle calculation.");
        }
        prev_TR = touching_right;
    }
    return std::min(prev_TR[N], prev_TR[Y]);  
}

int ViennaDangles::process_chains(
                    const std::vector<std::vector<size_t>>& dangle_chains,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies
                ) {
    int dangle_energy = 0;
    for (const std::vector<size_t>& chain : dangle_chains) {
        if (chain.size() == 1) {
            size_t idx = chain[0];
            int min_energy = dangle_energies[idx].min();
            dangle_energy += min_energy;
            continue;
        } else {
            dangle_energy += process_chain(chain, children, dangle_energies);
        }
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
    ML_Type ml_type; // Placeholder for future use
    bool front_dangle = children.front()->begin - node.begin <= 2;
    bool back_dangle  = node.end - children.back()->end <= 2;

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
    
    int dangle_energy = 0;
    return dangle_energy;
}

}