#pragma once
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>

#include "FileUtils.hpp"
#include "common.hpp"
using json = nlohmann::json;

namespace knotergy {

struct pk_param {
    pk_param()
        : name("DirksPierce09 pseudoknot params from HotKnotsV2"),
          pk_in_ext(-138),
          pk_in_mloop(1007),
          pk_in_pk(1500),
          band(246),
          unpaired_in_pk(6),
          cr_in_pk(96),
          pk_stack_x(0.89),
          pk_internal_x(0.74),
          pk_mloop_init(341),
          pk_mloop_bp(56),
          pk_mloop_unpaired(12) {}

    pk_param(const std::string& param_name, int pk_ext, int pk_multi, int pk_pk, int band_pen,
             int unpaired_pk, int cr_in_pk_pen, double pk_stack_multiplier,
             double pk_internal_multiplier, int pk_mloop_init_pen, int pk_mloop_bp_pen,
             int pk_mloop_unpaired_pen)
        : name(param_name),
          pk_in_ext(pk_ext),
          pk_in_mloop(pk_multi),
          pk_in_pk(pk_pk),
          band(band_pen),
          unpaired_in_pk(unpaired_pk),
          cr_in_pk(cr_in_pk_pen),
          pk_stack_x(pk_stack_multiplier),
          pk_internal_x(pk_internal_multiplier),
          pk_mloop_init(pk_mloop_init_pen),
          pk_mloop_bp(pk_mloop_bp_pen),
          pk_mloop_unpaired(pk_mloop_unpaired_pen) {}

    const std::string name;
    const int pk_in_ext;          // pseudoknot in exterior loop penalty
    const int pk_in_mloop;        // pseudoknot in multiloop penalty
    const int pk_in_pk;           // pseudoknot in pseudoloop penalty
    const int band;               // band penalty
    const int unpaired_in_pk;     // unpaired bases in pseudoknot penalty
    const int cr_in_pk;           // closed region nested in pseudoknot penalty
    const double pk_stack_x;      // stacked pair that spans a band penalty multiplier
    const double pk_internal_x;   // internal pair that spans a band penalty multiplier
    const int pk_mloop_init;      // multiloop that spans a band penalty
    const int pk_mloop_bp;        // base pair for multiloop that spans a band penalty
    const int pk_mloop_unpaired;  // unpaired bases in a multiloop that spans a band penalty
};

class PseudoknotParams {
   public:
    static inline std::shared_ptr<const pk_param> pkp = std::make_shared<const pk_param>();

    static const pk_param load_pk_param(
        const std::string& paramFile = "./params/common/pk_DirksPierce09_HotKnotsV2.json") {
        if (paramFile.empty()) {
            return {};
        }
        if (!FileUtils::file_exists(paramFile)) {
            THROW_ERROR("Pseudoknot parameters JSON file \"" + paramFile + "\" not found.");
        }

        if (!FileUtils::is_file(paramFile)) {
            THROW_ERROR("Pseudoknot parameters path \"" + paramFile + "\" is not a file.");
        }

        std::cout << "Loading pseudoknot parameters from: " << paramFile << std::endl;
        PseudoknotParams::pkp = std::make_shared<const pk_param>(parse_pk_json(paramFile));
        return *PseudoknotParams::pkp;
    }

    [[nodiscard]] static const pk_param parse_pk_json(const std::string& jsonFile) {
        std::ifstream f(jsonFile);
        if (!f.is_open()) {
            THROW_ERROR("Error: Unable to open pseudoknot parameter file: " + jsonFile);
        }
        json data = json::parse(f);
        json pk = data["pseudoknot_parameters"].get<json>();
        return pk_param(pk.value("name", "default"), pk.value("pk_in_ext", 0),
                        pk.value("pk_in_mloop", 0), pk.value("pk_in_pk", 0), pk.value("band", 0),
                        pk.value("unpaired_in_pk", 0), pk.value("cr_in_pk", 0),
                        pk.value("pk_stack_x", 1.0), pk.value("pk_internal_x", 1.0),
                        pk.value("pk_mloop_init", 0), pk.value("pk_mloop_bp", 0),
                        pk.value("pk_mloop_unpaired", 0));
    }
};

}  // namespace knotergy