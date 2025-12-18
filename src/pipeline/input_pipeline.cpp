#include "input_pipeline.hpp"


#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


#include "../energy/ComputeEnergy.hpp"
#include "../loop_tree/LoopFactory.hpp"
#include "../preprocessing/RNAEntry.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "shared.hpp"
#include "read_file.hpp"

namespace knotergy {


// Collects all RNA input entries from console and/or file.
std::vector<RNAEntry> get_all_inputs(const std::string& input_file, 
                                     const std::string& sequence,
                                     const std::string& structure) {
    std::vector<RNAEntry> entries;

    // get console input
    if (!sequence.empty()) {
        entries.emplace_back("Console Sequence", sequence, structure);
    }

    // get file inputs
    if (!input_file.empty()) {
        std::vector<RNAEntry> file_entries = get_all_file_entries(input_file);
        // move values into entries (avoids deep copies). Keeps console as first entry
        entries.insert(entries.end(), std::make_move_iterator(file_entries.begin()),
                       std::make_move_iterator(file_entries.end()));
    }
    if (entries.empty()) THROW_ERROR("No Input Data Given");
    return entries;
}

std::vector<ProcessedRNAEntry> process_inputs(const std::vector<RNAEntry>& inputs, const std::vector<modified_base_params>& modified_params) {
    std::vector<ProcessedRNAEntry> processed_inputs;
    processed_inputs.reserve(inputs.size());
    for (RNAEntry rna : inputs) {
        processed_inputs.emplace_back(RNAProcessor::process_rna(std::move(rna), modified_params));
    }
    return processed_inputs;
}

}  // namespace knotergy