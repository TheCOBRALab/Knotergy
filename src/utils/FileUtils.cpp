#include "FileUtils.hpp"

#include "utils/colors.hpp"

#include <iostream>
#include <stack>

namespace knotergy {

// Mostly used for resolving default data paths
std::string FileUtils::resolve_data_path(const std::string& path, bool check_exists) {
    // Conda package path
    const char* conda_prefix = std::getenv("CONDA_PREFIX");
    if (conda_prefix && *conda_prefix) {
        std::string conda_path =
            join_path(join_path(std::string(conda_prefix), "share/knotergy"), path);
        if (!check_exists || file_exists(conda_path)) {
            return conda_path;
        }
    }

    // Source-tree path (Knotergy is the root)
    std::string source_path = join_path(std::string(KNOTERGY_SOURCE_DIR), path);
    if (!check_exists || file_exists(source_path)) {
        return source_path;
    }

    // Explicit/direct path first (Likely won't be used, but just in case)
    if (!check_exists || file_exists(path)) {
        return path;
    }

    std::ostringstream msg;

    msg << "Data file not found: " << path << "\n"
        << "Checked source path: " << source_path << "\n"
        << "CONDA_PREFIX: " << (conda_prefix ? conda_prefix : "not set") << "\n"
        << "Checked conda path: "
        << (conda_prefix ? join_path(join_path(std::string(conda_prefix), "share/knotergy"), path)
                         : "not set");

    THROW_ERROR(msg.str());
}

// checks if the a slash should be added between base and path, and returns the joined path
std::string FileUtils::join_path(const std::string& base, const std::string& path) {
    if (base.empty()) return path;
    if (path.empty()) return base;

    if (base.back() == '/' || base.back() == '\\') {
        return base + path;
    }

    return base + "/" + path;
}

// filesystem::exists not supported in older macOS Conda packages
bool FileUtils::file_exists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

// Check if path is a file
bool FileUtils::is_file(const std::string& name) {
    struct stat buffer;
    if (stat(name.c_str(), &buffer) != 0) {
        return false;
    }
    return S_ISREG(buffer.st_mode);
}

// Check if path is a directory
bool FileUtils::is_directory(const std::string& name) {
    struct stat buffer;
    if (stat(name.c_str(), &buffer) != 0) {
        return false;
    }
    return S_ISDIR(buffer.st_mode);
}

// Get file modification time as uint64_t timestamp, or 0 if file doesn't exist or on error
std::uint64_t FileUtils::get_file_mtime(const std::string& path) {
    if (path.empty()) return 0;

    struct stat st = {};
    if (stat(path.c_str(), &st) != 0) return 0;

    return static_cast<std::uint64_t>(st.st_mtime);
}

// List files in a directory
std::vector<std::string> FileUtils::get_files_in_dir(const std::string& dir, bool include_dirs,
                                                     bool recursive) {
    std::vector<std::string> out;

    std::stack<std::string> dirs_to_process;
    dirs_to_process.push(dir);

    while (!dirs_to_process.empty()) {
        std::string current_dir = dirs_to_process.top();
        dirs_to_process.pop();

        DIR* d = ::opendir(current_dir.c_str());
        if (!d) {
            throw std::runtime_error("opendir failed: " + current_dir + " (" +
                                     std::strerror(errno) + ")");
        }

        // Read entries
        while (dirent* e = ::readdir(d)) {
            // skip . and ..
            if (std::strcmp(e->d_name, ".") == 0 || std::strcmp(e->d_name, "..") == 0) continue;

            // construct full path
            std::string full_path = current_dir;
            if (!full_path.empty() && full_path.back() != '/') full_path += '/';
            full_path += e->d_name;

            // add to list if it's a file
            if (is_file(full_path) || (include_dirs && is_directory(full_path)))
                out.push_back(full_path);
            if (recursive && is_directory(full_path)) dirs_to_process.push(full_path);
        }

        ::closedir(d);
    }

    return out;
}

// Helper to extract filename without extension from a path
std::string FileUtils::strip_extension(const std::string& path) {
    // 1. Extract filename
    size_t slash = path.find_last_of("/\\");
    std::string filename = (slash == std::string::npos) ? path : path.substr(slash + 1);

    // 2. Handle hidden files like ".bashrc"
    if (!filename.empty() && filename[0] == '.') {
        size_t dot = filename.find('.', 1);  // skip first dot
        if (dot == std::string::npos) {
            return filename;  // no extension
        }
    }

    // 3. Remove extension
    size_t dot = filename.find_last_of('.');
    if (dot != std::string::npos) {
        filename = filename.substr(0, dot);
    }

    return filename;
}

}  // namespace knotergy