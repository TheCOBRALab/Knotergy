#pragma once

#include <iomanip>
#include <iostream>

#include "PseudoknotParams.hpp"
#include "ViennaParams.hpp"

namespace knotergy {

class OutputManager {
   public:
    static void print_parameter_report(const vrna_md_param& vienna_params,
                                       const pk_param& pk_params, const all_mod_params& mp) {
        constexpr int W = 18;

        std::cout << "Parameters\n";
        std::cout << "------------------------------------------------\n";

        std::cout << std::left << std::setw(W)
                  << "ViennaRNA:" << vienna_params.get_source_info().resolved_name
                  << " (Dangles: " << vienna_params.md.dangles << ")\n";

        std::cout << std::left << std::setw(W)
                  << "Pseudoknot:" << pk_params.get_source_info().resolved_name << '\n';

        std::cout << std::left << std::setw(W)
                  << "Modified bases:" << std::to_string(mp.size()) + " loaded"
                  << (mp.empty() ? " (type -m to load)" : "") << '\n';

        // Only show details if explicitly requested AND not too many
        if (!mp.empty()) {
            std::cout << "\nModified bases\n";
            std::cout << "------------------------------------------------\n";

            for (const auto& param : mp.get_all_params()) {
                std::cout << "- " << param.name() << '\n';
            }
        }

        std::cout << '\n';
    }
};

}  // namespace knotergy