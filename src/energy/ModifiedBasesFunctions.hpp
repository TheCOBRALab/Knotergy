#pragma once

#include <utility>

#include "../preprocessing/RNAProcessor.hpp"
#include "./ViennaDangles.hpp"
#include "./ViennaFunctions.hpp"

namespace knotergy {

enum class ModLookup { Stacking, Terminal, Mismatch, Dangle5, Dangle3 };

// Stores the differences in energy contributions due to modified bases
struct ModDiffs {
   ModDiffs(int terminal_diff, int mismatch_diff, int n5d_diff, int n3d_diff):
         terminalAU{terminal_diff},
         mismatch{mismatch_diff},
         n5d{n5d_diff},
         n3d{n3d_diff} {}
    const int terminalAU;
    const int mismatch;
    const int n5d;
    const int n3d;
};

class ModifiedBasesFunctions {
   public:
    static int find_mod_stack_energy(const size_t& i, const size_t& j, const size_t& ci,
                                        const size_t& cj, std::string sequence,
                                        const std::vector<std::string_view>& mod_sequence,
                                        const std::vector<modified_base_params>& mod_params);
    
    static int find_mod_external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence,
                                           const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params);

   private:
    // Returns a vector of unique modified bases found at the specified indices
    static std::vector<std::string_view> unique_modified_bases_at_indices(std::vector<size_t> indices, const std::vector<std::string_view>& mod_sequence);

    // Joins string views at specified indices into a single string
    static std::string join_string_views(std::vector<size_t> indices,
                                         const std::vector<std::string_view>& mod_sequence);

    // Returns the modified energy based on the lookup type and key (returns unmodified energy if no modifications found)
    static int get_mod_energy(const std::string& key,
                                 const std::vector<std::string_view>& modified,
                                 const std::vector<modified_base_params>& mod_params,
                                 int unmod_energy, ModLookup lookup_type);
    
    static void modify_dangle_set( DangleSet& original_set, ModDiffs diffs);
   
   // Returns the difference between modified and unmodified energy based on the lookup type and key
    static int get_mod_energy_difference(const std::string& key,
                                          const std::vector<std::string_view>& modified,
                                          const std::vector<modified_base_params>& mod_params,
                                          int unmod_energy, ModLookup lookup_type);
   
   static ModDiffs get_mod_dangle_energy_diffs(const std::shared_ptr<LoopNode>& c, const int n5d,
                                               const int n3d, const unsigned int type, const unsigned int r_type,
                                               const std::vector<std::string_view>& unique_mod_bases,
                                               const std::vector<std::string_view>& mod_sequence,
                                               const std::vector<modified_base_params>& mod_params, bool is_external) ;
};

}  // namespace knotergy