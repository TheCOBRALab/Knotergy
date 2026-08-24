#pragma once

#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/ViennaParams.hpp"

#include <iomanip>
#include <iostream>

namespace knotergy {

class OutputManager {
   public:
    static constexpr std::string_view version() { return "0.3.1"; }

    static void print_banner() {
        std::cout << R"(                                    
  ▄▄▄▄   ▄▄▄                                      
 █▀ ██  ██               █▄                       
    ██ ██    ▄          ▄██▄      ▄        ▄      
    █████    ████▄ ▄███▄ ██ ▄█▀█▄ ████▄▄████ ██ ██
    ██ ██▄   ██ ██ ██ ██ ██ ██▄█▀ ██   ██ ██ ██▄██
  ▀██▀  ▀██▄▄██ ▀█▄▀███▀▄██▄▀█▄▄▄▄█▀   ▀████▄▄▀██▀
                                          ██   ██ 
                                        ▀▀▀  ▀▀▀  
)";
        std::cout << "Version " << version() << "\n\n";
    }

    static void print_version() { std::cout << "Knotergy " << version() << '\n'; }

    static void print_parameter_report(const vrna_md_param& vienna_params,
                                       const pk_param& pk_params, const all_mod_params& mp,
                                       bool efn2_correction) {
        constexpr int W = 18;

        print_version();

        std::cout << "------------------------------------------------\n";
        std::cout << "Parameters at " << vienna_params.md.temperature << "°C\n";

        std::cout << std::left << std::setw(W)
                  << "ViennaRNA:" << vienna_params.get_source_info().resolved_name
                  << " (Dangle Model: " << vienna_params.md.dangles << ")\n";

        std::cout << std::left << std::setw(W) << "Pseudoknot:" << pk_params.name
                  << " (Round Method: " << static_cast<int>(pk_params.round) << " ("
                  << pk_params.round << "))\n";

        std::cout << std::left << std::setw(W)
                  << "Modified bases:" << std::to_string(mp.size()) + " loaded"
                  << (mp.empty() ? " (use -m to load)" : "") << '\n';

        if (efn2_correction) {
            std::cout << std::left << std::setw(W) << "efn2 correction" << "Enabled\n";
        }

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