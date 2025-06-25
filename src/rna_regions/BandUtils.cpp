#include "BandUtils.hpp"

namespace compute_energy {

/*──────────────── extend a stacked helix by exactly one pair ─────────────*/
bool BandUtils::extend_stem(size_t& il, size_t& jr) {
    const size_t n = pairings_.size();
    size_t ip = il + 1;
    size_t jp = jr - 1;

    /* skip unpaired positions on either side */
    while (ip < n && pairings_[ip] == NULL_INDEX) ++ip;
    while (jp > 0 && pairings_[jp] == NULL_INDEX) --jp;

    if (ip < jp && pairings_[ip] == jp) {  // still a canonical stack
        il = ip;
        jr = jp;
        return true;
    }
    return false;
}

/*──────────────── find all bands inside [left, right] ────────────────────*/
std::vector<Band> BandUtils::find_bands_in_region(size_t left, size_t right) {
    std::vector<Band> bands;

    for (size_t i = left; i <= right; ++i) {
        /* skip anything that is   – already used
         *                         – unpaired
         *                         – closing half
         *                         – pairs outside this region */
        if (done_[i] || pairings_[i] == NULL_INDEX || pairings_[i] < i || pairings_[i] > right)
            continue;

        size_t j = pairings_[i];
        size_t il = i;
        size_t jr = j;

        // walks the stem until last base pair in the given region
        while (il > left && jr <= right && pairings_[il - 1] && extend_stem(il, jr)) {
        }

        bands.push_back(Band{i, il, jr, j});

        /* mark every position that belongs to this band */
        for (size_t k = i; k <= il; ++k) {
            done_[k] = true;
            if (pairings_[k] != NULL_INDEX) done_[pairings_[k]] = true;
        }
        for (size_t k = jr; k <= j; ++k) {
            done_[k] = true;
            if (pairings_[k] != NULL_INDEX) done_[pairings_[k]] = true;
        }

        i = il;  // fast-forward
    }
    return bands;
}

/*──────────────── attach bands + pseudo-nest info to one node ────────────*/
void BandUtils::annotate_bands(const std::shared_ptr<LoopNode>& node) {
    /* only pseudoknots need bands; leave the rest untouched */
    if (node->loop_type != LoopType::Pseudoknot) return;
    node->bands = find_bands_in_region(node->begin, node->end);
    node->number_of_bands = static_cast<int>(node->bands.size());

    /* classify every direct child */
    for (std::shared_ptr<LoopNode>& child : node->children) {
        child->pseudo_type = PseudoNestedType::OutsideBand;
        for (const Band& b : node->bands) {
            if (b.contains(child->begin) && b.contains(child->end)) {
                child->pseudo_type = PseudoNestedType::InsideBand;
                break;
            }
        }
    }
}
}  // namespace compute_energy