#pragma once
#include "../io/PseudoknotParams.hpp"
#include "../loop_tree/LoopNode.hpp"
#include "../preprocessing/ProcessedRNAEntry.hpp"
#include "../preprocessing/RNAProcessor.hpp"
#include "PseudoknotFunctions.hpp"
#include "ViennaFunctions.hpp"

namespace knotergy {
/**
 * @brief Computes the total energy of an RNA secondary structure including pseudoknots.
 *
 * This class processes a loop tree representation of an RNA structure and calculates
 * the Gibbs free energy contribution of each loop. It uses ViennaRNA for non-pseudoknotted
 * loops and specialized pseudoknot energy functions for pseudoknotted loops.
 */
class ComputeEnergy {
   public:
    /**
     * @brief Construct a ComputeEnergy object and process the loop tree.
     *
     * Upon construction, this immediately processes the entire loop tree starting
     * from the root node and calculates the total energy.
     *
     * @param root_node Shared pointer to the root node of the loop tree.
     * @param processed_rna The processed RNA entry containing sequence and structure information.
     * @param vienna_params ViennaRNA energy parameters.
     * @param pseudo_params Pseudoknot energy parameters.
     * @param mod_params Vector of modified base parameters for computing energy with modified
     * nucleotides.
     * @param round Whether to round energy values (default: false).
     * @param verbose Whether to print detailed energy breakdown (default: false).
     */
    ComputeEnergy(std::shared_ptr<LoopNode> root_node, const ProcessedRNAEntry& processed_rna,
                  vrna_md_param& vp, const knotergy::pk_param& pkp,
                  const std::vector<modified_base_param>& mp = {}, bool round = false,
                  bool verbose = false)
        : root_node_{root_node},
          processed_rna_{processed_rna},
          vp_{vp},
          pkp_{pkp},
          mp_{mp},
          sequence_{processed_rna.get_sequence()},
          round_{round} {
        process_tree(*root_node_, verbose);
    };

    /**
     * @brief Get the total computed energy of the RNA structure.
     *
     * @return The total Gibbs free energy in centicalories (hundredths of kcal/mol).
     */
    float getEnergy() const { return energy_; };
    bool getInfiniteEnergyFlag() const { return infinite_energy_flag_; };

   private:
    std::shared_ptr<LoopNode> root_node_;
    const ProcessedRNAEntry& processed_rna_;
    vrna_md_param& vp_;
    const knotergy::pk_param& pkp_;
    const std::vector<modified_base_param>& mp_;
    const std::string& sequence_;
    float energy_ = 0.0f;
    bool round_ = false;
    bool infinite_energy_flag_ =
        false;  ///< Flag to indicate if any loop has infinite energy (e.g., invalid structures).

    /**
     * @brief Process the entire loop tree and calculate energies.
     *
     * Recursively traverses the loop tree and computes energy for each node.
     *
     * @param root_node The root node of the loop tree to process.
     * @param verbose Whether to print verbose energy breakdown (default: false).
     */
    void process_tree(LoopNode& root_node, bool verbose = false);

    /**
     * @brief Process the entire loop tree with modified base energy calculations.
     *
     * Recursively traverses the loop tree and computes energy for each node,
     * using modified base parameters when applicable.
     *
     * @param root_node The root node of the loop tree to process.
     * @param verbose Whether to print verbose energy breakdown (default: false).
     */
    void process_modified_tree(LoopNode& node, bool verbose);

    /**
     * @brief Process a single loop node and compute its energy contribution.
     *
     * @param node The loop node to process.
     * @return The energy contribution of this node in centicalories.
     */
    float process_node(LoopNode& node);

    /**
     * @brief Process a single loop node with modified base energy calculations.
     *
     * This function computes the energy contribution of a loop node, taking into account
     * any modified bases present in the sequence and using the appropriate energy parameters.
     *
     * @param node The loop node to process.
     * @return The energy contribution of this node in centicalories.
     */
    float process_modified_node(LoopNode& node);
};

}  // namespace knotergy