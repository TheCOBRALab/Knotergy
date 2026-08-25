#pragma once
#include "energy/pseudoknots/PseudoknotFunctions.hpp"
#include "energy/vienna/ViennaFunctions.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
#include "preprocessing/RNAProcessor.hpp"

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
     * @param root_node Raw pointer to the root node of the loop tree.
     * @param processed_rna The processed RNA entry containing sequence and structure information.
     * @param vienna_params ViennaRNA energy parameters.
     * @param pseudo_params Pseudoknot energy parameters.
     * @param mod_params Vector of modified base parameters for computing energy with modified
     * nucleotides.
     * @param pk_dangles Whether to include pseudoknot dangle energy contributions (default: false).
     * @param efn2_correction Whether to apply the efn2 single-bulge correction (default: false).
     * @param verbose Whether to print detailed energy breakdown (default: false).
     */
    ComputeEnergy(LoopNode& root_node, const ProcessedRNAEntry& processed_rna, vrna_md_param& vp,
                  const knotergy::pk_param& pkp, const all_mod_params& mp = {},
                  bool pk_dangles = false, bool efn2_correction = false, bool verbose = false)
        : root_node_{root_node},
          pRNA_{processed_rna},
          vp_{vp},
          pkp_{pkp},
          mp_{mp},
          sequence_{processed_rna.get_sequence()},
          mod_sequence_{processed_rna.get_modified_sequence()},
          pk_dangles_{pk_dangles},
          efn2_correction_{efn2_correction},
          has_modified_bases_{processed_rna.has_modified_bases()} {
        process_tree(root_node_, verbose);
    };

    /**
     * @brief Get the total computed energy of the RNA structure.
     *
     * @return The total Gibbs free energy in centicalories (hundredths of kcal/mol).
     */
    [[nodiscard]] double getEnergy() const { return energy_; };
    [[nodiscard]] bool getInfiniteEnergyFlag() const { return infinite_energy_flag_; };

   private:
    LoopNode& root_node_;
    const ProcessedRNAEntry& pRNA_;
    vrna_md_param& vp_;
    const knotergy::pk_param& pkp_;
    const all_mod_params& mp_;
    const std::string& sequence_;
    const std::vector<std::string_view>& mod_sequence_;
    double energy_ = 0.0;
    bool infinite_energy_flag_ = false;  ///< True if any loop produces infinite energy
    bool pk_dangles_ = false;            ///< Include pseudoknot dangle energy contributions
    bool efn2_correction_ = false;       ///< Flag for efn2 single-bulge correction
    bool has_modified_bases_ = false;    ///< True if the RNA sequence contains modified bases

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
     * @brief Process a single loop node and compute its energy contribution.
     *
     * @param node The loop node to process.
     * @return The energy contribution of this node in centicalories.
     */
    double process_node(LoopNode& node);
};

}  // namespace knotergy