#pragma once

#include <array>
#include <iostream>
#include <stack>
#include <string>
#include <vector>
#include <unordered_map>

#include "../pipeline/shared.hpp"
#include "ClosedRegion.hpp"
#include "RNAEntry.hpp"
namespace knotergy {

class RNAProcessedEntry {
    public:
    explicit RNAProcessedEntry (
            RNAEntry rna,
            std::vector<size_t> pairings,
            std::vector<ClosedRegion> closed_regions,
            std::vector<size_t> closed_regions_pairings,
            std::vector<int> unpaired_prefix_sum
        ) : name_{rna.name},
            sequence_{rna.sequence},
            structure_{rna.structure},
            pairings_{pairings},
            closed_regions_{closed_regions},
            closed_regions_pairings_{closed_regions_pairings},
            unpaired_prefix_sum_{unpaired_prefix_sum} {}

            const std::string& get_name() const {return name_;}
            const std::string& get_sequence() const {return sequence_;}
            const std::string& get_structure() const {return structure_;}
            const std::vector<size_t>& get_pairings() const {return pairings_;}
            const std::vector<ClosedRegion>& get_closed_regions() const {return closed_regions_;}
            const std::vector<size_t>& get_closed_regions_pairings() const {return closed_regions_pairings_;}
            size_t size() const {return structure_.size();}

            int get_unpaired_count(size_t from, size_t to) const {
                if (from >= unpaired_prefix_sum_.size() || to > unpaired_prefix_sum_.size()) {
                    throw std::out_of_range("Index out of range in get_unpaired_count");
                }

                if (from >= to) return 0;

                // Return the difference in unpaired counts between the two indices
                return unpaired_prefix_sum_[to] - unpaired_prefix_sum_[from];
            }
            int get_unpaired_count(ClosedRegion closed_region) const {
                return get_unpaired_count(closed_region.begin, closed_region.end);
            }
    
    private : 
    const std::string name_;
    const std::string sequence_;
    const std::string structure_;
    const std::vector<size_t> pairings_;
    const std::vector<ClosedRegion> closed_regions_;
    const std::vector<size_t> closed_regions_pairings_;
    const std::vector<int> unpaired_prefix_sum_;

    

};

}  // namespace knotergy