#include "parser/parser_functions.hpp"

#include <dirent.h>
#include <sys/stat.h>

namespace {

std::string join_path(const std::string& dir, const std::string& name)
{
    if (dir.empty() || dir.back() == '/')
        return dir + name;
    return dir + '/' + name;
}

bool ends_with(const std::string& text, const std::string& suffix)
{
    return text.size() >= suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/* Adds every <name>.mol file under dir to found[name].  Symbolic links to
 * files count; symbolic links to directories are not followed, so a link
 * cycle cannot loop, and hidden directories (.git and the like) are skipped. */
void collect_mol_files(const std::string& dir, std::map<std::string, std::vector<std::string>>& found)
{
    DIR* dirStream { opendir(dir.c_str()) };
    if (dirStream == nullptr) {
        std::cerr << "Warning, cannot read directory " << dir << " while looking for mol files\n";
        return;
    }

    std::vector<std::string> subdirs {};
    while (const dirent* entry = readdir(dirStream)) {
        const std::string name { entry->d_name };
        if (name == "." || name == "..")
            continue;
        const std::string path { join_path(dir, name) };
        struct stat linkInfo {};
        if (lstat(path.c_str(), &linkInfo) != 0)
            continue;
        if (S_ISDIR(linkInfo.st_mode)) {
            if (name[0] != '.')
                subdirs.push_back(path);
        } else if (ends_with(name, ".mol") && name.size() > 4) {
            struct stat targetInfo {};
            if (stat(path.c_str(), &targetInfo) == 0 && S_ISREG(targetInfo.st_mode))
                found[name.substr(0, name.size() - 4)].push_back(path);
        }
    }
    closedir(dirStream);

    // readdir order depends on the file system; sort so messages are stable
    std::sort(subdirs.begin(), subdirs.end());
    for (const auto& subdir : subdirs)
        collect_mol_files(subdir, found);
}

[[noreturn]] void mol_file_paths_error(const std::string& inputFileName, int lineNumber, const std::string& message)
{
    std::cerr << "Error in the molFiles block of " << inputFileName << ", line " << lineNumber << ": " << message
              << ", exiting.\n";
    exit(1);
}

} // namespace

std::map<std::string, std::string> read_mol_file_paths(const std::string& inputFileName)
{
    std::map<std::string, std::string> molFilePaths {};
    std::ifstream inputFile { inputFileName };
    if (!inputFile)
        return molFilePaths; // the caller reports the unreadable file

    // relative paths are relative to the input file, not the working directory
    const auto slashPos = inputFileName.rfind('/');
    const std::string baseDir { slashPos == std::string::npos ? "" : inputFileName.substr(0, slashPos + 1) };

    bool inBlock { false };
    int blockStart { 0 };
    int lineNumber { 0 };
    std::string line;
    while (getline(inputFile, line)) {
        ++lineNumber;
        const std::string tmpLine { create_tmp_line(line) };
        if (!inBlock) {
            if (tmpLine == "startmolfiles") {
                inBlock = true;
                blockStart = lineNumber;
            }
            continue;
        }
        if (tmpLine == "endmolfiles") {
            inBlock = false;
            continue;
        }
        if (skipLine(tmpLine))
            continue;
        // another block's start line means this one was never closed
        if (tmpLine.compare(0, 5, "start") == 0 && tmpLine.find('=') == std::string::npos)
            mol_file_paths_error(inputFileName, blockStart, "start molFiles has no end molFiles");

        std::string entry { line };
        remove_comment(entry);
        const auto equalPos = entry.find('=');
        if (equalPos == std::string::npos)
            mol_file_paths_error(inputFileName, lineNumber, "expected molName = path");

        // keep only the characters the molecules block keeps in a name
        std::string molName {};
        for (const char c : entry.substr(0, equalPos))
            if (std::isalnum(static_cast<unsigned char>(c)))
                molName += c;

        std::string path { entry.substr(equalPos + 1) };
        const auto first = path.find_first_not_of(" \t\r");
        const auto last = path.find_last_not_of(" \t\r");
        path = first == std::string::npos ? "" : path.substr(first, last - first + 1);

        if (molName.empty() || path.empty())
            mol_file_paths_error(inputFileName, lineNumber, "expected molName = path");
        if (path[0] != '/')
            path = baseDir + path;
        if (!molFilePaths.emplace(molName, path).second)
            mol_file_paths_error(inputFileName, lineNumber, "a second path for " + molName);
    }
    if (inBlock)
        mol_file_paths_error(inputFileName, blockStart, "start molFiles has no end molFiles");

    return molFilePaths;
}

std::string locate_mol_file(const std::string& molName, const Parameters& params,
    const std::map<std::string, std::string>& molFilePaths)
{
    const auto mapped = molFilePaths.find(molName);
    if (mapped != molFilePaths.end())
        return mapped->second;

    const std::string fileName { molName + ".mol" };
    if (!params.molDirRecursive)
        return join_path(params.molDir, fileName);

    // One walk serves every molecule; the directory is fixed for the run.
    static std::string indexedDir {};
    static bool indexed { false };
    static std::map<std::string, std::vector<std::string>> found {};
    if (!indexed || indexedDir != params.molDir) {
        found.clear();
        collect_mol_files(params.molDir.empty() ? "." : params.molDir, found);
        indexedDir = params.molDir;
        indexed = true;
    }

    const std::string searchDir { params.molDir.empty() ? "." : params.molDir };
    const auto hit = found.find(molName);
    if (hit == found.end()) {
        std::cerr << "Cannot find " << fileName << " in " << searchDir << " or its subdirectories, exiting.\n";
        exit(1);
    }
    if (hit->second.size() > 1) {
        std::cerr << "Found " << fileName << " in more than one place under " << searchDir << ":\n";
        for (const auto& path : hit->second)
            std::cerr << "  " << path << '\n';
        std::cerr << "Name the one to use in a molFiles block of the input file, exiting.\n";
        exit(1);
    }
    return hit->second.front();
}

void warn_unused_mol_file_paths(const std::map<std::string, std::string>& molFilePaths,
    const std::vector<std::string>& providedMols, const std::string& inputFileName)
{
    for (const auto& entry : molFilePaths)
        if (std::find(providedMols.begin(), providedMols.end(), entry.first) == providedMols.end())
            std::cout << "Warning, the molFiles block of " << inputFileName << " names " << entry.first
                      << ", which is not in its molecules block\n";
}
