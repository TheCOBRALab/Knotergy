
#include "ViennaDangles.hpp"

namespace knotergy {

// State for if the previous pair in a chain took the right dangle
enum TouchingRight {
    RightFree = 0,
    RightTaken = 1
};

// Calculate dangle energies for external loops (dangle type 1) with precomputed dangle energies
int ViennaDangles::get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children, const std::vector<DangleSet>& dangle_energies) {
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    return process_chains(dangle_chains, children, dangle_energies);
}

// Calculate dangle energies for external loops (dangle type 1)
int ViennaDangles::get_external_dangle_1(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence) {
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, sequence);
    return get_external_dangle_1(children, dangle_energies);
}



// Calculate dangle energies for multibranch loops (dangle type 1)
int ViennaDangles::get_multi_dangle_1(const LoopNode& node, const std::string& sequence){
    const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
    bool is_external = false;
    std::vector<DangleSet> dangle_energies = populate_children_dangle_energies(children, sequence, is_external);
    DangleSet closing = get_ml_dangle_energy(node, sequence);
    std::vector<std::vector<size_t>> dangle_chains = get_dangle_chains(children);
    return process_ml_chains(dangle_chains, children, dangle_energies, node, closing);
}


DangleSet ViennaDangles::get_ml_dangle_energy(const LoopNode& node, const std::string& sequence){
    size_t pi = node.begin, pj = node.end;
    DangleSet ml_dangle; // closing pair dangles
    int n5d = vrna_nucleotide_encode(sequence[pi + 1], &ViennaParams::md);
    int n3d = vrna_nucleotide_encode(sequence[pj - 1], &ViennaParams::md);
    unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[pi], sequence[pj]);

    ml_dangle.no_dangle    = vrna_E_multibranch_stem(pair_type, -1,  -1,  ViennaParams::p);
    ml_dangle.left_dangle  = vrna_E_multibranch_stem(pair_type, -1, n5d,  ViennaParams::p);
    ml_dangle.right_dangle = vrna_E_multibranch_stem(pair_type, n3d, -1,  ViennaParams::p);
    ml_dangle.both_dangle  = vrna_E_multibranch_stem(pair_type, n3d, n5d, ViennaParams::p);
    return ml_dangle;
}

// Precompute dangle energies for all children in the loop
std::vector<DangleSet> ViennaDangles::populate_children_dangle_energies(
    const std::vector<std::shared_ptr<LoopNode>>& children,
    const std::string& sequence,
    const bool& is_external
) {
    // stores the 4 dangle energy options for each child
    std::vector<DangleSet> dangle_energies;

    // Function pointer to appropriate ViennaRNA dangle energy function
    int (*vrna_E_stem)(unsigned int, int, int, vrna_param_t*); 
    vrna_E_stem = is_external ? vrna_E_exterior_stem : vrna_E_multibranch_stem;
    
    // Compute dangle energies for each child
    for (size_t i = 0; i < children.size(); ++i) {
        
        // Get child node and its indices
        const std::shared_ptr<LoopNode>& child = children[i];
        const size_t& ci = child->begin, cj = child->end;

        // Convert child nucleotides to numerical encoding for ViennaRNA (-1 for out of bounds)
        int n5d = ci > 0 ? vrna_nucleotide_encode(sequence[ci - 1], &ViennaParams::md) : -1;
        int n3d = cj < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[cj + 1], &ViennaParams::md) : -1;
        unsigned int pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj]);

        // Store the four dangle energy options for this child
        DangleSet d_energy;
        d_energy.no_dangle    = vrna_E_stem(pair_type, -1,  -1,  ViennaParams::p);
        d_energy.left_dangle  = vrna_E_stem(pair_type, n5d, -1,  ViennaParams::p);
        d_energy.right_dangle = vrna_E_stem(pair_type, -1,  n3d, ViennaParams::p);
        d_energy.both_dangle  = vrna_E_stem(pair_type, n5d, n3d, ViennaParams::p);
        dangle_energies.push_back(d_energy);
    }
    return dangle_energies;
}

// Identify chains of children that share dangles 
std::vector<std::vector<size_t>> ViennaDangles::get_dangle_chains(const std::vector<std::shared_ptr<LoopNode>>& children) {
    std::vector<std::vector<size_t>> dangle_chains;

    // Stores child indices in each chain. Initialize with first child.
    dangle_chains.reserve(children.size());
    dangle_chains.push_back({0});

    // Iterate through children to identify chains
    for (size_t i = 1; i < children.size(); ++i) {
        const std::shared_ptr<LoopNode>& child = children[i];
        const std::shared_ptr<LoopNode>& prev_child = children[i - 1];
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }

        // Check if current child is contiguous with previous child (shares a dangle or is adjacent)
        if (child->begin - prev_child->end <= 2) {
            dangle_chains.back().push_back(i);
        } else{ // start new chain
            dangle_chains.push_back({i});
        }
    }
    return dangle_chains;
}

// Dynamic programming to compute optimal dangle energies for a single chain of children
int ViennaDangles::process_chain(
                    const std::vector<size_t>& chain,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const bool& disable_first_left_dangle,
                    const bool& disable_last_right_dangle,
                    std::array<int,2> init,
                    DangleSet closing
                ) {
    std::array<int,2> prev = init; // default {0, INF}
    for (size_t idx : chain) {
        const DangleSet& energies = dangle_energies[idx];
        std::array<int,2> cur = {INF, INF};

        // Check if left or right dangle is possible based on adjacency (no unpaired bases in between)
        bool disable_left_dangle = (idx != chain.front() && contiguous_children(*children[idx], *children[idx - 1])) || (idx == chain.front() && disable_first_left_dangle);
        bool disable_right_dangle = (idx != chain.back() && contiguous_children(*children[idx], *children[idx + 1])) || (idx == chain.back() && disable_last_right_dangle);
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
            cur[RightFree] = std::min({prev[RightFree] + eNone, prev[RightFree] + eLeft, prev[RightTaken] + eNone});
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
    return std::min(prev[RightFree] + closing.best(), prev[RightTaken] + closing.best_left());
}

// Process multiple chains of children and aggregate their dangle energies
int ViennaDangles::process_chains(
                    const std::vector<std::vector<size_t>>& dangle_chains,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const bool& disable_first_left_dangle,
                    const bool& disable_last_right_dangle,
                    std::array<int,2> init,
                    DangleSet closing
                ) {
    int dangle_energy = 0;
    size_t chain_idx = 0;
    for (const std::vector<size_t>& chain : dangle_chains) {
        if (dangle_chains.size() == 1) { // only one chain
            return process_chain(chain, children, dangle_energies, disable_first_left_dangle, disable_last_right_dangle, init, closing);
        } else if (chain_idx == 0) { // first chain
            dangle_energy += process_chain(chain, children, dangle_energies, disable_first_left_dangle, false, init);
        } else if (chain_idx == dangle_chains.size() - 1) { // last chain
            dangle_energy += process_chain(chain, children, dangle_energies, false, disable_last_right_dangle, {0, INF}, closing);
        } else { // middle chains
            dangle_energy += process_chain(chain, children, dangle_energies);
        }
        ++chain_idx;
    }
    return dangle_energy;
}

// Specialized processing for multibranch loops with closing pair dangles
int ViennaDangles::process_ml_chains(
                    const std::vector<std::vector<size_t>>& dangle_chains,
                    const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::vector<DangleSet>& dangle_energies,
                    const LoopNode& node,
                    const DangleSet closing
                ) {
    enum class ML_Type { None, Front, Back, Both, Loop};
   
    bool front_dangle = children.front()->begin - node.begin <= 2;
    bool back_dangle  = node.end - children.back()->end <= 2;
    std::cout << "ML Dangle Energies - No: " << closing.no_dangle 
              << ", Left: " << closing.left_dangle 
              << ", Right: " << closing.right_dangle 
              << ", Both: " << closing.both_dangle << std::endl;

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
        return closing.best() + process_chains(dangle_chains, children, dangle_energies);
    } 

    bool disable_first_left_dangle = false;
    bool disable_last_right_dangle = false;
    std::array<int,2> init = {0, INF};
    DangleSet closing_set = DangleSet();

    if (ml_type == ML_Type::Front){
        // ((....)..(....)...)
        // ^^ closing pair is touching first child
        if (contiguous(node.begin, children.front()->begin)){
            disable_first_left_dangle = true; // disables the left dangle on first child
            init = {closing.best_right(), INF}; // disables left dangle on closing pair
        } else {
            // (.(....)..(....)...)
            //  ^ closing pair exactly one unpaired base away from first child
            init = {closing.best_right(), closing.best()};
        }
    } else if (ml_type == ML_Type::Back){
        // (...(....)..((....))
        //                   ^^closing pair is touching last child
        if (contiguous(node.end, children.back()->end)){
            disable_last_right_dangle = true;
            init = {closing.best_left(), INF};
        } else {
            // (...(....)..(.....).)
            //                    ^ closing pair exactly one unpaired base away from last child
            closing_set = closing;
        }
    } else if (ml_type == ML_Type::Both || ml_type == ML_Type::Loop){
        // ((....)..(....))
        //  ^          ^
        if (contiguous(node.begin, children.front()->begin) && contiguous(node.end, children.back()->end)){
            disable_first_left_dangle = true;
            disable_last_right_dangle = true;
            init = {closing.no_dangle, INF};
        } else if (contiguous(node.begin, children.front()->begin)){
            // ((....)..(.....).)
            // ^^              ^ closing pair touching first child and last child is one base away
            disable_first_left_dangle = true;
            int chain1_energy = process_chains(
                                        dangle_chains,
                                        children,
                                        dangle_energies,
                                        disable_first_left_dangle,
                                        disable_last_right_dangle,
                                        {closing.no_dangle, INF}
                                    );
            disable_last_right_dangle = true;
            int chain2_energy = process_chains(
                                        dangle_chains,
                                        children,
                                        dangle_energies,
                                        disable_first_left_dangle,
                                        disable_last_right_dangle,
                                        {closing.best_right(), INF}
                                    );
            
            return std::min(chain1_energy, chain2_energy);
            
        } else if (contiguous(node.end, children.back()->end)){
            // (.(....)..(....))
            //  ^             ^^ closing pair one base away from first child and touching last child
            disable_last_right_dangle = true;
            init = {0, closing.best_left()};
        } else {
            // (.(....)..(.....).)
            //  ^               ^ closing pair one base away from both children
            int chain1_energy = process_chains(
                                        dangle_chains,
                                        children,
                                        dangle_energies,
                                        disable_first_left_dangle,
                                        disable_last_right_dangle,
                                        {0, closing.best_left()}
                                    );
            disable_last_right_dangle = true;
            int chain2_energy = process_chains(
                                        dangle_chains,
                                        children,
                                        dangle_energies,
                                        disable_first_left_dangle,
                                        disable_last_right_dangle,
                                        {0, closing.best()}
                                    );
            return std::min(chain1_energy, chain2_energy);
            
        }
    }

    return process_chains(
                dangle_chains,
                children,
                dangle_energies,
                disable_first_left_dangle,
                disable_last_right_dangle,
                init,
                closing_set
           );
    
    int dangle_energy = 0;
    return dangle_energy;
}
// ./build/Knotergy -p ./params/common/rna_turner2004.par -s AAAAAAAAAUUUUUUUUAAAAUUUUUUU -r "(((..(((....)))..((...)).)))" -d 1
// echo -e "AAAAAAAAAUUUUUUUUAAAAUUUUUUU\n(((..(((....)))..((...)).)))" | RNAeval -d 1

}

// echo -e "AAAAAAUUUUUUUAAAAAAUUUUUU\n((((....))....((....)).))" | RNAeval