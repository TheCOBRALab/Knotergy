#pragma once
#include "FileUtils.hpp"
#include "common.hpp"

#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace knotergy {

    struct pseudoknot_params {
        pseudoknot_params(const std::string& param_name = "", int pk_ext = 0, int pk_multi = 0, int pk_pk = 0,
                         int band_pen = 0, int unpaired_pk = 0, int nested_cr_pen = 0,
                         double pk_stack_multiplier = 1.0, double pk_internal_multiplier = 1.0,
                         int pk_mloop_init_pen = 0, int pk_mloop_bp_pen = 0, int pk_mloop_unpaired_pen = 0)
            : name(param_name),
              pk_in_ext(pk_ext),
              pk_in_mloop(pk_multi),
              pk_in_pk(pk_pk),
              band(band_pen),
              unpaired_in_pk(unpaired_pk),
              nested_cr(nested_cr_pen),
              pk_stack_x(pk_stack_multiplier),
              pk_internal_x(pk_internal_multiplier),
              pk_mloop_init(pk_mloop_init_pen),
              pk_mloop_bp(pk_mloop_bp_pen),
              pk_mloop_unpaired(pk_mloop_unpaired_pen)
        {}

        const std::string name;
        const int pk_in_ext;        // pseudoknot in exterior loop penalty
        const int pk_in_mloop;      // pseudoknot in multiloop penalty
        const int pk_in_pk;         // pseudoknot in pseudoloop penalty
        const int band;             // band penalty
        const int unpaired_in_pk;   // unpaired bases in pseudoknot penalty
        const int nested_cr;        // closed region nested in pseudoknot penalty
        const double pk_stack_x;    // stacked pair that spans a band penalty multiplier
        const double pk_internal_x; // internal pair that spans a band penalty multiplier
        const int pk_mloop_init;    // multiloop that spans a band penalty
        const int pk_mloop_bp;      // base pair for multiloop that spans a band penalty
        const int pk_mloop_unpaired; // unpaired bases in a multiloop that spans a band penalty
    };


    class PseudoknotParams {
        public:
        [[nodiscard]] static pseudoknot_params load_pk_params(const std::string& paramFile = ""){
                if (paramFile.empty()) {return {};}
                if (!FileUtils::file_exists(paramFile)) {
                    THROW_ERROR("Pseudoknot parameters JSON file \"" + paramFile + "\" not found.");
                }

                if (!FileUtils::is_file(paramFile)) {
                    THROW_ERROR("Pseudoknot parameters path \"" + paramFile + "\" is not a file.");
                }

                return parse_pk_json(paramFile);
        }

        [[nodiscard]] static pseudoknot_params parse_pk_json(const std::string& jsonFile) {
            std::ifstream f(jsonFile);
            if (!f.is_open()) {
                THROW_ERROR("Error: Unable to open pseudoknot parameter file: " + jsonFile);
            }
            json data = json::parse(f);
            json pk = data["pseudoknot_parameters"].get<json>();
            return pseudoknot_params(
                pk.value("name", ""),
                pk.value("pk_in_ext", 0),
                pk.value("pk_in_mloop", 0),
                pk.value("pk_in_pk", 0),
                pk.value("band", 0),
                pk.value("unpaired_in_pk", 0),
                pk.value("nested_cr", 0),
                pk.value("pk_stack_x", 1.0),
                pk.value("pk_internal_x", 1.0),
                pk.value("pk_mloop_init", 0),
                pk.value("pk_mloop_bp", 0),
                pk.value("pk_mloop_unpaired", 0)
            );
        }
    };

} // namespace knotergy