#pragma once

#include "../preprocessing/RNAProcessor.hpp"
#include "./ViennaFunctions.hpp"
#include "./ViennaDangles.hpp"
#include <utility>

namespace knotergy {

    class ModifiedBasesFunctions {
    public:
        static double find_mod_stack_energy(const size_t& i, const size_t& j, const size_t& ci, const size_t& cj, std::string sequence, const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params);

        static double find_mod_external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence, std::vector<std::string_view> mod_sequence, const std::vector<modified_base_params>& mod_params);

        static double find_mod_multibranch_energy(const LoopNode& node, const std::string& sequence, std::vector<std::string_view> mod_sequence, const std::vector<modified_base_params>& mod_params);
        
    private:
        // Returns a vector of unique modified bases found at the specified indices
        static std::vector<std::string_view> modified_bases(std::initializer_list<size_t> indices, const std::vector<std::string_view>& mod_sequence);

        // Joins string views at specified indices into a single string
        static std::string join_string_views(std::initializer_list<size_t> indices, const std::vector<std::string_view>& mod_sequence);

        // Returns the corrected energy for a given key and modified bases
        static double get_corrected_energy(const std::string& key, const std::vector<std::string_view>& modified, const std::vector<modified_base_params>& mod_params, double vienna_energy);


        static double get_corrected_branch_energy(const LoopNode& node, unsigned int type, int n5d, int n3d, std::vector<std::string_view> modified, const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params);


        static DangleSet correct_dangle_set(const DangleSet& original, const std::shared_ptr<LoopNode>& c, unsigned int type, int n5d, int n3d, std::vector<std::string_view> modified, const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params, bool is_external);
    };


} // namespace knotergy