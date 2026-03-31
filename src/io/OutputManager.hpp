#pragma once

#include <iomanip>
#include <iostream>

namespace knotergy {

class OutputManager {
   public:
    static void print_parameter_report(const vrna_md_param& vienna_params,
                                       const pk_param& pk_params,
                                       const std::vector<modified_base_param>& modified_params) {
        constexpr int W = 18;

        std::cout << "Parameters\n";
        std::cout << "------------------------------------------------\n";

        std::cout << std::left << std::setw(W)
                  << "ViennaRNA:" << vienna_params.get_source_info().resolved_name << '\n';

        std::cout << std::left << std::setw(W)
                  << "Pseudoknot:" << pk_params.get_source_info().resolved_name << '\n';

        std::cout << std::left << std::setw(W)
                  << "Modified bases:" << std::to_string(modified_params.size()) + " loaded"
                  << (modified_params.empty() ? " (type -m to load)" : "") << '\n';

        // Only show details if explicitly requested AND not too many
        if (!modified_params.empty()) {
            std::cout << "\nModified bases\n";
            std::cout << "------------------------------------------------\n";

            for (const auto& param : modified_params) {
                std::cout << "- " << param.name << '\n';
            }
        }

        std::cout << '\n';
    }
};

}  // namespace knotergy