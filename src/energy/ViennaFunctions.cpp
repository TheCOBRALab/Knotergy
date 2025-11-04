
#include "ViennaFunctions.hpp"
#include <algorithm>

namespace knotergy {

const vrna_param_t* ViennaFunctions::get_parameters() const {return P;};

int ViennaFunctions::stack_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence){
    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        std::cerr << "Invalid indices for stack energy calculation." << std::endl;
        return 0;
    }
    // c = child or nested base pair
    unsigned int type1 = get_pair_type(sequence[i], sequence[j]);
    unsigned int type2 = get_pair_type(sequence[ci], sequence[cj]);
    return P->stack[type1][reverse_pair_type(type2)];
}

int ViennaFunctions::stack_energy(BasePair pair, BasePair child, const std::string& sequence) {
    return stack_energy(pair.i, pair.j, child.i, child.j, sequence);
}

int ViennaFunctions::hairpin_energy(size_t i, size_t j, const std::string& sequence) {
    if (j <= i || j >= sequence.size()) {
        std::cerr << "Invalid indices for hairpin energy calculation." << std::endl;
        return 0;
    }

    // loop size
    unsigned int size = static_cast<unsigned int>(j - i - 1);
    if (size < 3) {
        std::cerr << "Warning: Hairpin loop size is less than 3. Infinite Energy. Sequence: "
                  << sequence << " i: " << i << ", j: " << j << std::endl;
    }

    unsigned int pair_type = get_pair_type(sequence[i], sequence[j]);
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &md);
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &md);

    // If loop size < 7, you MUST pass the loop sequence substring
    // https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/eval/hairpin.h#L78C1-L81C94
    const char* loop_seq = nullptr;
    std::string loop_subseq;
    if (size < 7) {
        loop_subseq = sequence.substr(i, size + 2);  // length == size
        loop_seq = loop_subseq.c_str();              // temporary c-string
    }

    return vrna_E_hairpin(size, pair_type, si1, sj1, loop_seq, P);
}

int ViennaFunctions::hairpin_energy(const BasePair& pair, const std::string& sequence) {
    return hairpin_energy(pair.i, pair.j, sequence);
}

int ViennaFunctions::internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence) {
    // c = child or nested bp

    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        THROW_ERROR("Invalid indices for internal loop energy calculation."); // Add i, j, ci, cj to message
    }

    unsigned int n1 = static_cast<unsigned int>(ci - i - 1);
    unsigned int n2 = static_cast<unsigned int>(j - cj - 1);
    unsigned int pair_type = get_pair_type(sequence[i], sequence[j]);
    unsigned int pair_type2 = reverse_pair_type(get_pair_type(sequence[ci], sequence[cj]));
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &md);   // 5' mismatch of closing pair
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &md);   // 3' mismatch of closing pair
    int sp1 = vrna_nucleotide_encode(sequence[ci - 1], &md);  // 5' mismatch of enclosed pair
    int sq1 = vrna_nucleotide_encode(sequence[cj + 1], &md);  // 3' mismatch of enclosed pair

    return vrna_E_internal(n1, n2, pair_type, pair_type2, si1, sj1, sp1, sq1, P);
}

int ViennaFunctions::internal_loop_energy(BasePair pair, BasePair child, const std::string& sequence) {
    return internal_loop_energy(pair.i, pair.j, child.i, child.j, sequence);
}

int ViennaFunctions::multibranch_energy(const LoopNode& node, const std::string& sequence) {
    size_t i = node.begin;
    size_t j = node.end;

    // penalties
    int energy = P->MLclosing; // closing penalty
    energy += node.exclusive_unpaired_bases_count * P->MLbase; // unpaired bases penalty

    int n5d, n3d;
    switch (md.dangles){
        case 0:
            n5d = -1;
            n3d = -1;
            break;
        case 1:
            // return get_multi_dangle_1(node, sequence) + energy;
        case 2:
            n5d = vrna_nucleotide_encode(sequence[i + 1], &md);
            n3d = vrna_nucleotide_encode(sequence[j - 1], &md);
            break;
        default:
            THROW_ERROR("Dangle model `" + std::to_string(md.dangles) + "` not supported in multi-branch energy calculation.");
    }

    unsigned int pair_type = get_pair_type(sequence[i], sequence[j]);

    // energy of the closing pair
    energy += vrna_E_multibranch_stem(reverse_pair_type(pair_type), n3d, n5d, P);

    const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
    for (const std::shared_ptr<LoopNode>& child : children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            // energy += child->number_of_bands * P->MLintern[child_pair_type];
            continue;
        }

        size_t ci = child->begin;
        size_t cj = child->end;
        unsigned int child_pair_type = get_pair_type(sequence[ci], sequence[cj]);

        int child_n5d;
        int child_n3d;
        switch (md.dangles){
            case 0:
                child_n5d = -1;
                child_n3d = -1;
                break;
            case 2:
                child_n5d = vrna_nucleotide_encode(sequence[ci - 1], &md);
                child_n3d = vrna_nucleotide_encode(sequence[cj + 1], &md);
                break;
            default:
                // Dangle 1 is handled separately above
                THROW_ERROR("Dangle model `" + std::to_string(md.dangles) + "` not supported in multi-branch energy calculation. (HOW DID WE GET HERE?)");
        }
        energy += vrna_E_multibranch_stem(child_pair_type, child_n5d, child_n3d, P);
    }

    return energy;
}

int ViennaFunctions::external_energy(const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::string& sequence) {

    if (md.dangles == 1) {
        return get_external_dangle_1(children, sequence);
    }

    int energy = 0;
    for (std::shared_ptr<LoopNode> c : children) {
        if (c->loop_type != LoopType::Pseudoknot) {
            unsigned int pair_type = get_pair_type(sequence[c->begin], sequence[c->end]);
            int n5d = c->begin > 0 ? vrna_nucleotide_encode(sequence[c->begin - 1], &md) : -1;
            int n3d = c->end < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[c->end + 1], &md) : -1;
            if (md.dangles == 0) {
                n5d = -1;
                n3d = -1;
            }
            energy += vrna_E_exterior_stem(pair_type, n5d, n3d, P);
        }
    }
    return energy;
}

unsigned int ViennaFunctions::get_pair_type(const char& i, const char& j) {
    int encoded_i = vrna_nucleotide_encode(i, &md);
    int encoded_j = vrna_nucleotide_encode(j, &md);
    return vrna_get_ptype_md(encoded_i, encoded_j, &md);
}

// I created a pull request to remove the need for static cast
// https://github.com/ViennaRNA/ViennaRNA/pull/270
unsigned int ViennaFunctions::reverse_pair_type(unsigned int type) const {
    return static_cast<unsigned int>(md.rtype[type]);
}

int ViennaFunctions::get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence) {
    std::vector<std::array<int,4>> dangle_energies = populate_children_dangle_energies(children, sequence);
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);

    enum TouchingRight {N = 0, Y = 1}; // No, Yes
    int dangle_energy = 0;
    // choose best dangle configuration for each child
    for (const std::vector<size_t>& chain : dangle_chains) {
        if (chain.size() == 1) {
            size_t idx = chain[0];
            int min_energy = *std::min_element(dangle_energies[idx].begin(), dangle_energies[idx].end());
            dangle_energy += min_energy;
            continue;
        }  else {
            // Dynamic programming over the chain to find optimal dangle configuration
            // prev_TR[0] = no right dangle on previous (Not touching right)
            // prev_TR[1] = right dangle on previous (Touching right))
            // touching_right[0] = no right dangle on current
            // touching_right[1] = right dangle on current
            std::array<int,2> prev_TR = {0, 0};
            for (size_t idx : chain) {
                const std::array<int,4>& energies = dangle_energies[idx];
                std::array<int,2> touching_right = {INF, INF};

                // Check if left or right dangle is possible based on adjacency (no unpaired bases in between)
                bool no_left_dangle = idx != chain.front() && (children[idx]->begin - children[idx - 1]->end == 1);
                bool no_right_dangle = idx != chain.back() && (children[idx + 1]->begin - children[idx]->end == 1);

                // energies: no dangle, left dangle, right dangle, both dangles
                int eNone = energies[None], eLeft = energies[Left], eRight = energies[Right], eBoth = energies[Both];

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
            dangle_energy += std::min(prev_TR[N], prev_TR[Y]);    
        }
    }
    return dangle_energy;
}


int ViennaFunctions::get_multi_dangle_1(const LoopNode& node, const std::string& sequence) {
    if (node.children.empty()) {
        THROW_ERROR("Multi-branch loop has no children for dangle calculation. (WTF?)");
    }

    const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
    std::array<int,4> ml_dangle = populate_ml_dangle_energies(node, sequence);
    std::vector<std::array<int,4>> dangle_energies = populate_children_dangle_energies(children, sequence, false);
    
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    // std::array<int,2> prev_TR = {std::min(ml_dangle[None], ml_dangle[Left]), std::min(ml_dangle[Right], ml_dangle[Both])};
}

std::array<int,4> ViennaFunctions::populate_ml_dangle_energies(const LoopNode& node, const std::string& sequence) {
    size_t pi = node.begin, pj = node.end;
    std::array<int,4> ml_dangle; // closing pair dangles
    int n5d = vrna_nucleotide_encode(sequence[pi + 1], &md);
    int n3d = vrna_nucleotide_encode(sequence[pj - 1], &md);
    unsigned int pair_type = get_pair_type(sequence[pi], sequence[pj]);
    unsigned int rev_pair_type = reverse_pair_type(pair_type);

    ml_dangle[None]  = vrna_E_multibranch_stem(rev_pair_type, -1,  -1,  P);
    ml_dangle[Left]  = vrna_E_multibranch_stem(rev_pair_type, -1,  n5d, P); 
    ml_dangle[Right] = vrna_E_multibranch_stem(rev_pair_type, n3d, -1,  P); 
    ml_dangle[Both]  = vrna_E_multibranch_stem(rev_pair_type, n3d, n5d, P);
    return ml_dangle;
}

std::vector<std::array<int,4>> ViennaFunctions::populate_children_dangle_energies(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence, const bool& is_external) {
    std::vector<std::array<int,4>> dangle_energies;
    int (*vrna_E_stem)(unsigned int, int, int, vrna_param_t*); 
    vrna_E_stem = is_external ? vrna_E_exterior_stem : vrna_E_multibranch_stem;

    for (size_t i = 0; i < children.size(); ++i) {
        const std::shared_ptr<LoopNode>& child = children[i];
        const size_t& ci = child->begin, cj = child->end;

        // Compute all dangle energies for this child
        int n5d = ci > 0 ? vrna_nucleotide_encode(sequence[ci - 1], &md) : -1;
        int n3d = cj < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[cj + 1], &md) : -1;
        unsigned int pair_type = get_pair_type(sequence[ci], sequence[cj]);

        std::array<int,4> d_energy;
        d_energy[None]  = vrna_E_stem(pair_type, -1,  -1,  P);
        d_energy[Left]  = vrna_E_stem(pair_type, n5d, -1,  P);
        d_energy[Right] = vrna_E_stem(pair_type, -1,  n3d, P);
        d_energy[Both]  = vrna_E_stem(pair_type, n5d, n3d, P);
        dangle_energies.push_back(d_energy);
    }
    return dangle_energies;
}

std::vector<std::vector<size_t>> ViennaFunctions::get_dangle_chains(const std::vector<std::shared_ptr<LoopNode>>& children) {
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


} // namespace knotergy