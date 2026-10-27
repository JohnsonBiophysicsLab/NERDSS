#include "classes/class_Parameters.hpp"
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <sys/stat.h>

namespace {

void print_usage(const char *program, std::ostream &os) {
  os << "Usage: " << program << " -f <parmfile> [options]\n"
     << "       " << program << " -r <restartfile> [options]\n"
     << "\n"
     << "  -f, --parmfile <file>    parameter file\n"
     << "  -s, --seed <n>           RNG seed, a whole number from 0 to "
     << UINT_MAX << "\n"
     << "  -r, --restart <file>     restart from this file\n"
     << "  -a, --add <file>         add file applied on restart\n"
     << "  -c, --coordinate <file>  initial coordinates\n"
     << "  -m, --mol-dir <dir>      read <name>.mol files from this directory\n"
     << "  -M, --mol-dir-recursive <dir>\n"
     << "                           ... or from anywhere under it\n"
     << "  -v, --verbose            more log output; -vv for the most\n"
     << "      --debug-force-dissoc force dissociation whenever possible\n"
     << "      --debug-force-assoc  force association whenever possible\n"
     << "      --print-system-info  print the system information\n"
     << "  -h, --help               show this message and exit\n";
}

// strtoul alone accepts "-1" (wrapping it), trailing junk and values above
// UINT_MAX on LP64, so check all three.
bool parse_seed(const char *text, unsigned &seed) {
  if (text[0] == '\0' || text[0] == '-' || text[0] == '+')
    return false;
  errno = 0;
  char *end{nullptr};
  unsigned long value{std::strtoul(text, &end, 10)};
  if (errno == ERANGE || *end != '\0' || value > UINT_MAX)
    return false;
  seed = static_cast<unsigned>(value);
  return true;
}

bool is_directory(const char *path) {
  struct stat info {};
  return stat(path, &info) == 0 && S_ISDIR(info.st_mode);
}

} // namespace

void parse_command(int argc, char *argv[], Parameters &params,
                   std::string &paramFileName, std::string &restartFileName,
                   std::string &addFileName, std::string &coordinateFileName,
                   unsigned int &seed) {
  std::cout << "Command:";
  for (int argItr{0}; argItr < argc; ++argItr)
    std::cout << ' ' << argv[argItr];
  std::cout << '\n' << std::flush;

  // Values for the long-only options, beyond the range of any short option
  enum : int {
    OPT_DEBUG_FORCE_DISSOC = 256,
    OPT_DEBUG_FORCE_ASSOC,
    OPT_PRINT_SYSTEM_INFO,
  };
  static const struct option longOptions[] = {
      {"parmfile", required_argument, nullptr, 'f'},
      {"seed", required_argument, nullptr, 's'},
      {"restart", required_argument, nullptr, 'r'},
      {"add", required_argument, nullptr, 'a'},
      {"coordinate", required_argument, nullptr, 'c'},
      {"verbose", no_argument, nullptr, 'v'},
      {"help", no_argument, nullptr, 'h'},
      {"debug-force-dissoc", no_argument, nullptr, OPT_DEBUG_FORCE_DISSOC},
      {"debug-force-assoc", no_argument, nullptr, OPT_DEBUG_FORCE_ASSOC},
      {"print-system-info", no_argument, nullptr, OPT_PRINT_SYSTEM_INFO},
      {"mol-dir", required_argument, nullptr, 'm'},
      {"mol-dir-recursive", required_argument, nullptr, 'M'},
      {nullptr, 0, nullptr, 0}};

  bool molDirGiven{false};
  int opt{};
  while ((opt = getopt_long(argc, argv, "f:s:r:a:c:m:M:vh", longOptions,
                            nullptr)) != -1) {
    switch (opt) {
    case 'f':
      paramFileName = optarg;
      break;
    case 's':
      if (!parse_seed(optarg, seed)) {
        std::cerr << "Error reading seed '" << optarg << "', exiting.\n";
        exit(1);
      }
      break;
    case 'r':
      // nerdss_mpi used to open restart.dat<rank> whatever the argument, by
      // testing Parameters::rank, which nothing initializes, so which file
      // either build opened depended on what that int happened to hold.  A
      // restart file is one file, written by either build under the same
      // name (nerdss_mpi's by write_gathered_restart()); read the one named.
      restartFileName = optarg;
      params.fromRestart = true;
      break;
    case 'a':
      addFileName = optarg;
      break;
    case 'c':
      coordinateFileName = optarg;
      break;
    case 'v': // -v is 1, -vv (or -v -v) is 2
      if (params.debugParams.verbosity < 2)
        ++params.debugParams.verbosity;
      break;
    case 'h':
      print_usage(argv[0], std::cout);
      exit(0);
    case OPT_DEBUG_FORCE_DISSOC:
      params.debugParams.forceDissoc = true;
      break;
    case OPT_DEBUG_FORCE_ASSOC:
      params.debugParams.forceAssoc = true;
      break;
    case OPT_PRINT_SYSTEM_INFO:
      params.debugParams.printSystemInfo = true;
      break;
    case 'm':
    case 'M': {
      const bool recursive{opt == 'M'};
      if (molDirGiven && recursive != params.molDirRecursive) {
        std::cerr << "-m/--mol-dir and -M/--mol-dir-recursive cannot be "
                     "combined, exiting.\n";
        exit(1);
      }
      if (!is_directory(optarg)) {
        std::cerr << "Mol file directory '" << optarg
                  << "' is not a directory, exiting.\n";
        exit(1);
      }
      params.molDir = optarg;
      params.molDirRecursive = recursive;
      molDirGiven = true;
      break;
    }
    default: // getopt_long has already said what was wrong
      print_usage(argv[0], std::cerr);
      exit(1);
    }
  }

  if (optind < argc) {
    std::cerr << "Unexpected argument '" << argv[optind] << "', exiting.\n";
    print_usage(argv[0], std::cerr);
    exit(1);
  }
}
