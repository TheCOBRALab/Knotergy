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

/**
 * @brief Utility functions for file system operations.
 *
 * Provides cross-platform file system operations without relying on C++17 filesystem
 * (which is not fully supported in older macOS Conda packages).
 */
class FileUtils {
   public:
    /**
     * @brief Check if a file or directory exists.
     *
     * @param name Path to check.
     * @return True if the path exists.
     */
    [[nodiscard]] static bool file_exists(const std::string& name);

    /**
     * @brief Check if a path points to a regular file.
     *
     * @param name Path to check.
     * @return True if the path is a regular file.
     */
    [[nodiscard]] static bool is_file(const std::string& name);

    /**
     * @brief Check if a path points to a directory.
     *
     * @param name Path to check.
     * @return True if the path is a directory.
     */
    [[nodiscard]] static bool is_directory(const std::string& name);

    /**
     * @brief Get a list of files in a directory.
     *
     * @param dir Directory path to list.
     * @param include_dirs Whether to include subdirectories in the result (default: false).
     * @param recursive Whether to recursively list files in subdirectories (default: false).
     * @return Vector of file paths.
     */
    [[nodiscard]] static std::vector<std::string> get_files_in_dir(const std::string& dir,
                                                                   bool include_dirs = false,
                                                                   bool recursive = false);
};

}  // namespace knotergy