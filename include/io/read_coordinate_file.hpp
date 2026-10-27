/*! \file read_coordinate_file.hpp
 * \ingroup IO
 * \brief Reads the initial-coordinate file given with `-c` / `--coordinate`.
 *
 * Three formats are accepted.  Each data row names one molecule and gives its
 * center of mass in nm, in NERDSS's internal frame (the box runs from -L/2 to
 * +L/2 on every axis; NERDSS's own PDB and XYZ outputs are shifted by +L/2, so
 * an output frame has to be shifted back before it is fed in):
 *
 *  - PDB: `ATOM` records whose atom name (columns 13-16) is `COM`; the
 *    molecule name is columns 18-21 and x/y/z are columns 31-54.  Every other
 *    line is ignored.  A PDB residue name holds at most 4 characters.
 *  - XYZ: line 1 is the number of rows, line 2 is a comment and is always
 *    skipped (even when empty), then `<name> <x> <y> <z>` rows separated by
 *    whitespace.  Fewer rows than line 1 announces is an error; extra rows
 *    are ignored with a warning.
 *  - CSV: `<name>,<x>,<y>,<z>` rows, with an optional header.  The first
 *    non-blank row is a header when its 2nd or 3rd field is not entirely a
 *    number.  Spaces around fields, CRLF line endings and blank lines are
 *    tolerated.
 *
 * The format is taken from the extension (`.pdb`, `.xyz`, `.csv`, case
 * folded).  Any other or missing extension is sniffed: a line starting with
 * `ATOM` or `HETATM` means PDB, a first non-blank line that is a single
 * integer means XYZ, and a first non-blank line holding a comma means CSV.
 *
 * A malformed XYZ or CSV row stops the run with the file name and the line
 * number.  The PDB reader keeps its original behaviour and does not validate.
 */
#pragma once

#include "classes/class_Vec3D.hpp"

#include <map>
#include <string>
#include <vector>

enum class CoordinateFileFormat { pdb, xyz, csv };

/*! \ingroup IO
 * \brief The format's name as printed in the log ("PDB", "XYZ", "CSV").
 */
const char* coordinate_file_format_name(CoordinateFileFormat format);

/*! \ingroup IO
 * \brief Picks the format of a coordinate file from its extension, or from
 * its content when the extension decides nothing.  Exits with a message when
 * neither does.  `how` (optional) receives "extension" or "content".
 */
CoordinateFileFormat detect_coordinate_file_format(const std::string& fileName, std::string* how = nullptr);

/*! \ingroup IO
 * \brief Reads a coordinate file in the given format.
 * \return name -> COM positions, in file order.  Names are the molecule
 * names exactly as written; nothing is matched against the system here.
 */
std::map<std::string, std::vector<Vec3D>> read_coordinate_file(const std::string& fileName,
    CoordinateFileFormat format);

/*! \ingroup IO
 * \brief Detects the format, prints it, and reads the file.
 */
std::map<std::string, std::vector<Vec3D>> read_coordinate_file(const std::string& fileName);
