#include <dirent.h>  // For directory iteration in POSIX systems
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "io/io.hpp"

struct ComplexKey {
  std::string components;
  bool operator==(const ComplexKey& other) const {
    return components == other.components;
  }
};

namespace std {
template <>
struct hash<ComplexKey> {
  std::size_t operator()(const ComplexKey& k) const {
    return std::hash<std::string>()(k.components);
  }
};
}  // namespace std

void processFile(const std::string& filename,
                 std::map<double, std::unordered_map<ComplexKey, int>>& data) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return;
  }

  std::string line;
  double currentTime = -1;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    if (line.find("Time (s):") != std::string::npos) {
      // New time block
      std::string dummy;
      iss >> dummy >> dummy >> currentTime;
    } else {
      // Complex data
      int count;
      iss >> count;
      std::string components, component;
      while (iss >> component) {
        components += component + " ";
      }
      data[currentTime][{components}] += count;
    }
  }
  file.close();
}


// ---- the per-rank tallies, merged under the serial build's names ----
//
// Every rank writes its own DATA/<name>_<rank>.dat on the serial schedule
// (src/io_mpi/write_output.cpp), as its own share: bound pairs, event
// counters, observables and the transition matrix are counts the rank's own
// reactions made, the monomer and dimer counts are of the complexes it owns,
// the association log is the reactions it performed.  So the system's file
// is the sum of the ranks' files row by row (observables: minus nprocs - 1
// times the start value every rank was given), the transition matrix's
// lifetime lists laid end to end, and the association log sorted by step.
// Every rank writes the same rows on the same steps, so the files are merged
// line by line; a row missing on one rank ends the merge of that file with a
// note on stderr rather than a misaligned sum.
namespace {

bool open_rank_files(const std::string& prefix, int n,
                     std::vector<std::ifstream>& files) {
  files.clear();
  files.resize(n);
  for (int i = 0; i < n; ++i) {
    files[i].open(prefix + std::to_string(i) + ".dat");
    if (!files[i].is_open()) {
      std::cerr << "merge_outputs: no " << prefix << i << ".dat" << std::endl;
      return false;
    }
  }
  return true;
}

// Reads one line from every rank; false at the first end of file.  A rank
// whose file ends while others go on is reported once.
bool read_row(std::vector<std::ifstream>& files, std::vector<std::string>& rows,
              const std::string& what) {
  rows.assign(files.size(), "");
  int got = 0;
  for (size_t i = 0; i < files.size(); ++i)
    if (std::getline(files[i], rows[i])) ++got;
  if (got == static_cast<int>(files.size())) return true;
  if (got != 0)
    std::cerr << "merge_outputs: the ranks' " << what
              << " files differ in length; merged up to the shortest" << std::endl;
  return false;
}

// Fields of a separated row, empties kept (a row may end in a tab).
std::vector<std::string> split_fields(const std::string& line, char sep) {
  std::vector<std::string> fields;
  std::string field;
  std::istringstream iss(line);
  while (std::getline(iss, field, sep)) fields.push_back(field);
  if (!line.empty() && line.back() == sep) fields.push_back("");
  return fields;
}

std::string join_fields(const std::vector<std::string>& fields, char sep) {
  std::string line;
  for (size_t i = 0; i < fields.size(); ++i) {
    if (i) line += sep;
    line += fields[i];
  }
  return line;
}

bool is_integer(const std::string& s) {
  if (s.empty()) return false;
  size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
  if (i == s.size()) return false;
  for (; i < s.size(); ++i)
    if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
  return true;
}

// A header line copied from rank 0, then rows whose first field (the time)
// is rank 0's and whose other fields are summed over the ranks.  Integers
// stay integers; a field that is not a number (an empty trailing field) is
// copied from rank 0.
void merge_summed_rows(const std::string& prefix, const std::string& outName,
                       int n, char sep) {
  std::vector<std::ifstream> files;
  if (!open_rank_files(prefix, n, files)) return;
  std::ofstream out(outName);
  if (!out.is_open()) {
    std::cerr << "merge_outputs: cannot write " << outName << std::endl;
    return;
  }
  std::vector<std::string> rows;
  if (!read_row(files, rows, prefix)) return;
  out << rows[0] << '\n';
  while (read_row(files, rows, prefix)) {
    std::vector<std::string> merged = split_fields(rows[0], sep);
    std::vector<long long> sums(merged.size(), 0);
    std::vector<double> dsums(merged.size(), 0.0);
    std::vector<char> numeric(merged.size(), 1), integral(merged.size(), 1);
    for (int r = 0; r < n; ++r) {
      std::vector<std::string> fields = split_fields(rows[r], sep);
      if (r && !fields.empty() && !merged.empty() && fields[0] != merged[0])
        std::cerr << "merge_outputs: " << prefix << r << ".dat has row " << fields[0]
                  << " where rank 0 has " << merged[0] << std::endl;
      for (size_t j = 1; j < merged.size() && j < fields.size(); ++j) {
        if (is_integer(fields[j])) {
          sums[j] += std::stoll(fields[j]);
          dsums[j] += std::stod(fields[j]);
        } else {
          char* end = nullptr;
          const double v = std::strtod(fields[j].c_str(), &end);
          if (end && *end == '\0' && !fields[j].empty()) {
            integral[j] = 0;
            dsums[j] += v;
          } else {
            numeric[j] = 0;
          }
        }
      }
    }
    for (size_t j = 1; j < merged.size(); ++j) {
      if (!numeric[j]) continue;
      if (integral[j]) {
        merged[j] = std::to_string(sums[j]);
      } else {
        std::ostringstream oss;
        oss << dsums[j];
        merged[j] = oss.str();
      }
    }
    out << join_fields(merged, sep) << '\n';
  }
}

// The association event counters: blocks headed "time (s): <t>", then one
// "<label>: <count>" line per non-zero entry, in the serial writer's order
// (3D, then 3D to 2D, then 2D; dimers first, then by size).  Summed per
// label across the ranks; a label any rank wrote is written.
struct EventLabelOrder {
  bool operator()(const std::string& a, const std::string& b) const {
    return key(a) < key(b);
  }
  static std::pair<int, long long> key(const std::string& label) {
    int category = 3;
    if (label.compare(0, 6, "3Dto2D") == 0)
      category = 1;
    else if (label.compare(0, 2, "3D") == 0)
      category = 0;
    else if (label.compare(0, 2, "2D") == 0)
      category = 2;
    long long size = 0;  // "dimers"
    const size_t pos = label.find('n');
    if (pos != std::string::npos && pos + 1 < label.size() &&
        (label[pos + 1] == '=' || label[pos + 1] == '>')) {
      size = std::atoll(label.c_str() + pos + 2);
      // "n>K" is the open-ended last entry, after "n=K to ..."
      if (label[pos + 1] == '>') size = size * 2 + 1;
    }
    return {category, size};
  }
};

void merge_event_counters(const std::string& prefix, const std::string& outName, int n) {
  std::vector<std::ifstream> files;
  if (!open_rank_files(prefix, n, files)) return;
  std::ofstream out(outName);
  if (!out.is_open()) return;
  // Every rank writes the same block headers in the same order; within a
  // block the labels differ, so each rank is read a block at a time.
  std::vector<std::string> pendingHeader(n);
  bool primed = false;
  while (true) {
    std::string header;
    std::map<std::string, long long, EventLabelOrder> block;
    int withHeader = 0;
    for (int r = 0; r < n; ++r) {
      std::string line;
      std::string head = primed ? pendingHeader[r] : "";
      if (!primed) {
        while (std::getline(files[r], line))
          if (line.compare(0, 9, "time (s):") == 0) {
            head = line;
            break;
          }
      }
      pendingHeader[r].clear();
      if (head.empty()) continue;
      ++withHeader;
      if (header.empty()) header = head;
      else if (head != header)
        std::cerr << "merge_outputs: " << prefix << r << ".dat has block '" << head
                  << "' where rank 0 has '" << header << "'" << std::endl;
      while (std::getline(files[r], line)) {
        if (line.compare(0, 9, "time (s):") == 0) {
          pendingHeader[r] = line;  // the next block's header, kept for the next round
          break;
        }
        const size_t colon = line.rfind(": ");
        if (colon == std::string::npos) continue;
        block[line.substr(0, colon)] += std::atoll(line.c_str() + colon + 2);
      }
    }
    primed = true;
    if (withHeader == 0) break;
    if (withHeader != n) {
      std::cerr << "merge_outputs: the ranks' " << prefix
                << " files differ in their blocks; merged up to the shortest" << std::endl;
      break;
    }
    out << header << '\n';
    for (const auto& entry : block) out << entry.first << ": " << entry.second << '\n';
  }
}

// The observables: a header from rank 0, rows "t,v1,v2,...".  Every rank's
// observablesList starts from the whole system's value and counts its own
// reactions on top, so the system's value is the sum less nprocs - 1 start
// values (as write_gathered_restart() computes it).  A row at time 0 is
// the parser's literal "0,0" and is copied from rank 0.
void merge_observables(const std::string& prefix, const std::string& outName, int n,
                       const std::map<std::string, int>& observablesAtStart) {
  std::vector<std::ifstream> files;
  if (!open_rank_files(prefix, n, files)) return;
  std::ofstream out(outName);
  if (!out.is_open()) return;
  // The start values in the writer's order (write_observables() walks the
  // map).  A new run's file opens with a header line; a restart's does not
  // (the serial build appends to the run's file), so the header is optional.
  std::vector<long long> start;
  for (const auto& observable : observablesAtStart) start.push_back(observable.second);
  std::vector<std::string> rows;
  bool first = true;
  while (read_row(files, rows, prefix)) {
    if (first) {
      first = false;
      if (rows[0].compare(0, 4, "Time") == 0) {
        out << rows[0] << '\n';
        continue;
      }
    }
    std::vector<std::string> merged = split_fields(rows[0], ',');
    if (merged.empty()) continue;
    if (merged[0] == "0") {
      out << rows[0] << '\n';
      continue;
    }
    std::vector<long long> sums(merged.size(), 0);
    for (int r = 0; r < n; ++r) {
      std::vector<std::string> fields = split_fields(rows[r], ',');
      for (size_t j = 1; j < merged.size() && j < fields.size(); ++j)
        sums[j] += std::atoll(fields[j].c_str());
    }
    for (size_t j = 1; j < merged.size(); ++j)
      merged[j] = std::to_string(
          sums[j] - (j - 1 < start.size() ? (n - 1) * start[j - 1] : 0));
    out << join_fields(merged, ',') << '\n';
  }
}

// The transition matrix file: blocks of "time: <t>", the matrix of every
// counted template under its name (rows of integers, summed), then the
// lifetime lists under "lifetime for each mol type:", one "size of the
// cluster:<k>" line and one row of values per size (laid end to end in rank
// order).  Every rank writes the same lines in the same order, only the
// numbers differ, so the files are merged line by line.
void merge_transition_matrix(const std::string& prefix, const std::string& outName, int n) {
  std::vector<std::ifstream> files;
  if (!open_rank_files(prefix, n, files)) return;
  std::ofstream out(outName);
  if (!out.is_open()) return;
  std::vector<std::string> rows;
  bool lifetimes = false;
  while (read_row(files, rows, prefix)) {
    const std::string& line = rows[0];
    if (line.compare(0, 5, "time:") == 0) {
      lifetimes = false;
      out << line << '\n';
      continue;
    }
    if (line.compare(0, 8, "lifetime") == 0) {
      lifetimes = true;
      out << line << '\n';
      continue;
    }
    const bool valueRow = line.empty() || line[0] == ' ';
    if (!valueRow) {  // a section header, a template name, "size of the cluster:k"
      out << line << '\n';
      continue;
    }
    if (lifetimes) {
      for (int r = 0; r < n; ++r) out << rows[r];
      out << '\n';
      continue;
    }
    std::vector<long long> sums;
    for (int r = 0; r < n; ++r) {
      std::istringstream iss(rows[r]);
      long long v;
      size_t j = 0;
      while (iss >> v) {
        if (j >= sums.size()) sums.push_back(0);
        sums[j++] += v;
      }
    }
    for (long long v : sums) out << ' ' << v;
    out << '\n';
  }
}

// The association/dissociation log ("ITR:<step>,BOND|BREAK,..."), every
// rank's lines in step order, rank order within a step.
void merge_assoc_dissoc(const std::string& prefix, const std::string& outName, int n) {
  std::vector<std::ifstream> files;
  if (!open_rank_files(prefix, n, files)) return;
  std::ofstream out(outName);
  if (!out.is_open()) return;
  std::vector<std::pair<long long, std::string>> lines;
  for (int r = 0; r < n; ++r) {
    std::string line;
    while (std::getline(files[r], line)) {
      if (line.empty()) continue;
      long long step = 0;
      if (line.compare(0, 4, "ITR:") == 0) step = std::atoll(line.c_str() + 4);
      lines.emplace_back(step, line);
    }
  }
  std::stable_sort(lines.begin(), lines.end(),
                   [](const std::pair<long long, std::string>& a,
                      const std::pair<long long, std::string>& b) { return a.first < b.first; });
  for (const auto& entry : lines) out << entry.second << '\n';
}

}  // namespace

void merge_outputs(int totalrank, int molTemplateNum, bool mergePdbFiles,
                   const std::map<std::string, int>& observablesAtStart) {
  // MERGE THE PDB FILES IN THE PDB FOLDER AND STORE THE MERGED FILES IN
  // MERGEPDB FOLDER
  int n = totalrank;
  int skiplines = 23;
  std::string pdb_folder = "PDB";
  std::string merge_folder = "mergePDB";

  // COLLECT ALL .PDB FILENAMES
  std::vector<std::string> pdb_filenames;

  DIR* dir;
  struct dirent* ent;
  if ((dir = opendir(pdb_folder.c_str())) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      std::string filename = ent->d_name;
      if (filename.find(".cif") != std::string::npos) {
        pdb_filenames.push_back(filename);
      }
    }
    closedir(dir);
  } else if (mergePdbFiles) {
    // pdbWrite <= 0 makes no PDB folder, so its absence is only an error
    // when PDB files were written
    std::cerr << "Could not open directory: " << pdb_folder << std::endl;
  }

  // EXTRACT UNIQUE TIMEFRAMES
  std::set<int> timeframes;
  for (const auto& filename : pdb_filenames) {
    size_t underscorePos = filename.find("_");
    size_t dotPos = filename.find(".");
    if (underscorePos != std::string::npos && dotPos != std::string::npos) {
      int timeframe = std::stoi(
          filename.substr(underscorePos + 1, dotPos - underscorePos - 1));
      timeframes.insert(timeframe);
    }
  }

  // MERGE
  for (int timeframe : timeframes) {
    std::string merged_filename = std::to_string(timeframe) + ".cif";
    std::string merged_path = merge_folder + "/" + merged_filename;

    std::ofstream merged_file(merged_path);
    if (!merged_file.is_open()) {
      std::cerr << "Failed to open " << merged_path << std::endl;
      continue;
    }

    for (int rank = 0; rank < n; ++rank) {
      std::string input_filename =
          std::to_string(rank) + "_" + std::to_string(timeframe) + ".cif";
      std::string input_path = pdb_folder + "/" + input_filename;

      std::ifstream input_file(input_path);
      if (!input_file.is_open()) {
        std::cerr << "Failed to open " << input_path << std::endl;
        continue;
      }

      if (rank == 0) {
        merged_file << input_file.rdbuf();
      } else {
        std::string line;
        for (int i = 0; i < skiplines; ++i) {
          std::getline(input_file, line);  // Skip lines
        }
        merged_file << input_file.rdbuf();
      }
      input_file.close();
    }
    merged_file.close();
  }

  // MERGE THE COPY NUMBER FILES FROM ALL THE PROCESSORS AND STORE IN THE
  // MERGEOUT FOLDER
  std::string input_prefix = "DATA/copy_numbers_time_";
  std::string output_folder = "mergeOUT";
  std::string output_filename = output_folder + "/copy_numbers_time.dat";
  // Open output file
  std::ofstream output_file(output_filename);
  if (!output_file.is_open()) {
    std::cerr << "Failed to open output file for writing." << std::endl;
  }

  // PROCESS EACH FILE
  std::string line;
  std::vector<std::ifstream> input_files(n);
  bool headerProcessed = false;
  int numColumns = 0;
  for (int i = 0; i < n; ++i) {
    input_files[i].open(input_prefix + std::to_string(i) + ".dat");
    if (!input_files[i].is_open()) {
      std::cerr << "Failed to open input file " << i << std::endl;
    }

    if (!headerProcessed) {
      if (std::getline(input_files[i], line)) {
        output_file << line << std::endl;
        numColumns = std::count(line.begin(), line.end(), ',');
        headerProcessed = true;
      }
    } else {
      std::getline(input_files[i], line);
    }
  }

  while (true) {
    std::map<double, std::vector<double>> summed_values;
    bool end_of_files = false;
    for (int i = 0; i < n; ++i) {
      if (std::getline(input_files[i], line)) {
        std::istringstream iss(line);
        std::string value;
        double time;
        getline(iss, value, ',');
        time = std::stod(value);
        if (summed_values.find(time) == summed_values.end()) {
          summed_values[time] = std::vector<double>(numColumns, 0.0);
        }
        for (int j = 0; j < numColumns; ++j) {
          if (getline(iss, value, ',')) {
            summed_values[time][j] += std::stod(value);
          }
        }
      } else {
        end_of_files = true;
        break;
      }
    }
    if (end_of_files) break;
    for (const auto& entry : summed_values) {
      output_file << entry.first;
      for (const auto& sum : entry.second) {
        output_file << "," << sum;
      }
      output_file << std::endl;
    }
  }
  for (auto& file : input_files) {
    file.close();
  }
  output_file.close();

  // MERGE HISTOGRAM FILES FROM ALL PROCESSORS
  input_prefix = "DATA/histogram_complexes_time_";
  output_filename = output_folder + "/histogram_complexes_time.dat";

  std::map<double, std::unordered_map<ComplexKey, int>> mergedData;

  for (int i = 0; i < n; ++i) {
    processFile(input_prefix + std::to_string(i) + ".dat", mergedData);
  }

  std::ofstream output_file_hist(output_filename);
  if (!output_file_hist.is_open()) {
    std::cerr << "Failed to open output file for writing." << std::endl;
  }

  for (const auto& timeData : mergedData) {
    output_file_hist << "Time (s): " << timeData.first << std::endl;
    for (const auto& complexCount : timeData.second) {
      output_file_hist << complexCount.second << "\t"
                       << complexCount.first.components << std::endl;
    }
  }
  output_file_hist.close();

  // THE OTHER PER-RANK TALLIES, UNDER THE SERIAL BUILD'S NAMES
  merge_summed_rows("DATA/bound_pair_time_", output_folder + "/bound_pair_time.dat", n, '\t');
  merge_summed_rows("DATA/mono_dimer_time_", output_folder + "/mono_dimer_time.dat", n, '\t');
  merge_event_counters("DATA/event_counters_time_", output_folder + "/event_counters_time.dat", n);
  merge_observables("DATA/observables_time_", output_folder + "/observables_time.dat", n,
                    observablesAtStart);
  merge_transition_matrix("DATA/transition_matrix_time_",
                          output_folder + "/transition_matrix_time.dat", n);
  merge_assoc_dissoc("DATA/assoc_dissoc_time_", output_folder + "/assoc_dissoc_time.dat", n);
}
