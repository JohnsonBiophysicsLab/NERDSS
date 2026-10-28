/*! \file read_coordinate_file.cpp
 * \ingroup IO
 * \brief Readers for the `-c` coordinate file: PDB, XYZ and CSV.
 *
 * See read_coordinate_file.hpp for the formats and the detection rules.
 */
#include "io/read_coordinate_file.hpp"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {

using CoordinateMap = std::map<std::string, std::vector<Vec3D>>;

std::string trim(const std::string& s)
{
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return {};
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

std::string lower(std::string s)
{
    for (char& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

/* True when the whole of `s` (after trimming) is one floating-point number.
 * strtod alone accepts a numeric prefix ("1.5abc" -> 1.5), which is what a
 * header cell like "x" or "1st" must not pass as. */
bool parse_number(const std::string& s, double& value)
{
    const std::string t = trim(s);
    if (t.empty())
        return false;
    char* end = nullptr;
    value = std::strtod(t.c_str(), &end);
    return end == t.c_str() + t.size();
}

/* True when `s` is a non-negative base-10 integer with nothing else on it. */
bool parse_count(const std::string& s, long& value)
{
    const std::string t = trim(s);
    if (t.empty())
        return false;
    for (char c : t)
        if (!std::isdigit(static_cast<unsigned char>(c)))
            return false;
    value = std::strtol(t.c_str(), nullptr, 10);
    return true;
}

[[noreturn]] void fail_row(const std::string& fileName, long lineNo, const std::string& what, const std::string& line)
{
    std::cerr << "Error: coordinate file " << fileName << ", line " << lineNo << ": " << what << "\n"
              << "  " << line << "\n";
    exit(1);
}

/* Stores `name x y z` given as four already-split fields, or fails the run. */
void store_row(CoordinateMap& out, const std::vector<std::string>& fields, const std::string& fileName,
    long lineNo, const std::string& line)
{
    if (fields.size() != 4) {
        fail_row(fileName, lineNo,
            "expected 4 fields (name x y z), found " + std::to_string(fields.size()), line);
    }
    const std::string name = trim(fields[0]);
    if (name.empty())
        fail_row(fileName, lineNo, "empty molecule name", line);
    double xyz[3];
    const char* axis[3] = { "x", "y", "z" };
    for (int k = 0; k < 3; ++k) {
        if (!parse_number(fields[k + 1], xyz[k]))
            fail_row(fileName, lineNo, std::string(axis[k]) + " coordinate '" + trim(fields[k + 1]) + "' is not a number", line);
    }
    out[name].push_back(Vec3D { xyz[0], xyz[1], xyz[2] });
}

std::vector<std::string> split_whitespace(const std::string& line)
{
    std::vector<std::string> fields;
    std::istringstream iss(line);
    std::string f;
    while (iss >> f)
        fields.push_back(f);
    return fields;
}

std::vector<std::string> split_commas(const std::string& line)
{
    std::vector<std::string> fields;
    std::string::size_type start = 0;
    for (;;) {
        const auto comma = line.find(',', start);
        fields.push_back(trim(line.substr(start, comma == std::string::npos ? std::string::npos : comma - start)));
        if (comma == std::string::npos)
            return fields;
        start = comma + 1;
    }
}

/* The original -c reader, unchanged: ATOM records whose atom name is COM. */
CoordinateMap read_pdb(std::ifstream& inputFile)
{
    CoordinateMap moleculeCoordinates;
    std::string line;
    while (std::getline(inputFile, line)) {
        // Skip lines that don't start with "ATOM"
        if (line.substr(0, 4) != "ATOM") {
            continue;
        }

        // Parse the line according to PDB format
        std::string atomType = line.substr(12, 4);
        // Trim whitespace from atomType
        atomType.erase(0, atomType.find_first_not_of(" "));
        atomType.erase(atomType.find_last_not_of(" ") + 1);

        // Only process COM atom type (center of mass)
        if (atomType != "COM") {
            continue;
        }

        std::string moleculeName = line.substr(17, 4);
        // Trim whitespace from moleculeName
        moleculeName.erase(0, moleculeName.find_first_not_of(" "));
        moleculeName.erase(moleculeName.find_last_not_of(" ") + 1);

        // Extract x, y, z coordinates
        double x = std::stod(line.substr(30, 8));
        double y = std::stod(line.substr(38, 8));
        double z = std::stod(line.substr(46, 8));

        // Store coordinates for this molecule type in a vector
        // Each entry in the vector represents one molecule of that type from the PDB
        moleculeCoordinates[moleculeName].push_back(Vec3D { x, y, z });
    }
    return moleculeCoordinates;
}

CoordinateMap read_xyz(std::ifstream& inputFile, const std::string& fileName)
{
    CoordinateMap out;
    std::string line;
    long lineNo = 0;

    if (!std::getline(inputFile, line))
        fail_row(fileName, 1, "empty file; line 1 must be the number of rows", "");
    ++lineNo;
    long declared = 0;
    if (!parse_count(line, declared))
        fail_row(fileName, lineNo, "line 1 must be the number of rows, a whole number", line);

    // Line 2 is the comment and is skipped whatever it holds, even nothing.
    if (!std::getline(inputFile, line)) {
        if (declared == 0)
            return out;
        fail_row(fileName, 2, "file ends before the comment line", "");
    }
    ++lineNo;

    long rows = 0;
    while (std::getline(inputFile, line)) {
        ++lineNo;
        const std::vector<std::string> fields = split_whitespace(line);
        if (fields.empty())
            continue; // blank line
        if (rows == declared) {
            std::cout << "WARNING: coordinate file " << fileName << " announces " << declared
                      << " rows on line 1 but holds more; the rows from line " << lineNo
                      << " on are ignored\n";
            break;
        }
        store_row(out, fields, fileName, lineNo, line);
        ++rows;
    }
    if (rows < declared) {
        std::cerr << "Error: coordinate file " << fileName << " announces " << declared
                  << " rows on line 1 but holds only " << rows << "\n";
        exit(1);
    }
    return out;
}

CoordinateMap read_csv(std::ifstream& inputFile, const std::string& fileName)
{
    CoordinateMap out;
    std::string line;
    long lineNo = 0;
    bool first = true;
    while (std::getline(inputFile, line)) {
        ++lineNo;
        if (trim(line).empty())
            continue;
        const std::vector<std::string> fields = split_commas(line);
        if (first) {
            first = false;
            // A header is a first row whose 2nd or 3rd field is not a number.
            double v;
            const bool numeric = fields.size() >= 3 && parse_number(fields[1], v) && parse_number(fields[2], v);
            if (!numeric)
                continue;
        }
        store_row(out, fields, fileName, lineNo, line);
    }
    return out;
}

} // namespace

const char* coordinate_file_format_name(CoordinateFileFormat format)
{
    switch (format) {
    case CoordinateFileFormat::pdb:
        return "PDB";
    case CoordinateFileFormat::xyz:
        return "XYZ";
    case CoordinateFileFormat::csv:
        return "CSV";
    }
    return "?";
}

CoordinateFileFormat detect_coordinate_file_format(const std::string& fileName, std::string* how)
{
    const auto slash = fileName.find_last_of("/\\");
    const std::string base = slash == std::string::npos ? fileName : fileName.substr(slash + 1);
    const auto dot = base.find_last_of('.');
    if (dot != std::string::npos && dot > 0) {
        const std::string ext = lower(base.substr(dot + 1));
        if (how)
            *how = "extension";
        if (ext == "pdb")
            return CoordinateFileFormat::pdb;
        if (ext == "xyz")
            return CoordinateFileFormat::xyz;
        if (ext == "csv")
            return CoordinateFileFormat::csv;
    }

    if (how)
        *how = "content";
    std::ifstream inputFile(fileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open coordinate file: " << fileName << std::endl;
        exit(1);
    }
    // ATOM/HETATM anywhere wins (a PDB usually opens with HEADER, REMARK or
    // CRYST1); otherwise the first non-blank line decides.
    std::string line;
    std::string firstNonBlank;
    while (std::getline(inputFile, line)) {
        const std::string t = trim(line);
        if (t.compare(0, 4, "ATOM") == 0 || t.compare(0, 6, "HETATM") == 0)
            return CoordinateFileFormat::pdb;
        if (firstNonBlank.empty() && !t.empty())
            firstNonBlank = t;
    }
    long n;
    if (parse_count(firstNonBlank, n))
        return CoordinateFileFormat::xyz;
    if (firstNonBlank.find(',') != std::string::npos)
        return CoordinateFileFormat::csv;

    std::cerr << "Error: cannot tell the format of coordinate file " << fileName
              << ": the extension is not .pdb, .xyz or .csv, no line is an ATOM/HETATM record, "
                 "the first line is not a row count and no row is comma-separated\n";
    exit(1);
}

std::map<std::string, std::vector<Vec3D>> read_coordinate_file(const std::string& fileName,
    CoordinateFileFormat format)
{
    std::ifstream inputFile(fileName);
    if (!inputFile.is_open()) {
        // The original reader warned and carried on with the generated
        // coordinates; that is kept for every format.
        std::cerr << "Error: Could not open coordinate file: " << fileName << std::endl;
        return {};
    }
    switch (format) {
    case CoordinateFileFormat::pdb:
        return read_pdb(inputFile);
    case CoordinateFileFormat::xyz:
        return read_xyz(inputFile, fileName);
    case CoordinateFileFormat::csv:
        return read_csv(inputFile, fileName);
    }
    return {};
}

std::map<std::string, std::vector<Vec3D>> read_coordinate_file(const std::string& fileName)
{
    if (!std::ifstream(fileName).is_open()) {
        std::cerr << "Error: Could not open coordinate file: " << fileName << std::endl;
        return {};
    }
    std::string how;
    const CoordinateFileFormat format = detect_coordinate_file_format(fileName, &how);
    std::cout << "coordinate file format: " << coordinate_file_format_name(format) << " (from " << how << ")"
              << std::endl;
    return read_coordinate_file(fileName, format);
}
