#pragma once

#include <dirent.h>
#include <sys/stat.h>

#include <fstream>
#include <stack>
#include <string>
#include <vector>

#include "../preprocessing/RNAEntry.hpp"
#include "common.hpp"

namespace knotergy {

class FileUtils {
   public:
    // filesystem::exists not supported in older macOS Conda packages
    [[nodiscard]] static bool file_exists(const std::string& name);

    // Check if path is a file
    [[nodiscard]] static bool is_file(const std::string& name);

    // Check if path is a directory
    [[nodiscard]] static bool is_directory(const std::string& name);

    // Get files in a directory
    [[nodiscard]] static std::vector<std::string> get_files_in_dir(const std::string& dir,
                                                                   bool include_dirs = false,
                                                                   bool recursive = false);
};

}  // namespace knotergy