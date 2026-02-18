<p align="center">
  <img src="docs/nerdss_banner_export.svg" alt="Project banner" width="100%">
</p>


## NERDSS Development

<p align="center">
  <!-- Release -->
  <a href="https://github.com/JohnsonBiophysicsLab/NERDSS/releases">
    <img alt="Release" src="https://img.shields.io/github/v/release/JohnsonBiophysicsLab/NERDSS">
  </a>
  <!-- License -->
  <a href="https://github.com/JohnsonBiophysicsLab/NERDSS/blob/master/LICENSE">
    <img alt="License" src="https://img.shields.io/github/license/JohnsonBiophysicsLab/NERDSS">
  </a>
  <!-- C++ standard -->
  <img alt="C++" src="https://img.shields.io/badge/C%2B%2B-%E2%89%A514-blue">
  <!-- CMake minimum -->
  <img alt="Make" src="https://img.shields.io/badge/Make-%E2%89%A53.79-informational">
  <!-- Package managers (optional) -->
  <!-- img alt="Conan" src="https://img.shields.io/badge/Conan-ready-0ea5e9">
  <!-- img alt="vcpkg" src="https://img.shields.io/badge/vcpkg-port-22c55e">
  <!-- Platforms -->
  <img alt="Platforms" src="https://img.shields.io/badge/Linux%20|%20macOS%20-supported-success">
  <!-- Code style -->
  <!-- img alt="clang-format" src="https://img.shields.io/badge/clang--format-enforced-brightgreen"-->
</p>

Structure-Resolved Reaction-Diffusion Simulation Software by Johnson Lab, JHU

### Installation

To build serial NERDSS, you need:

1. A C++ compiler:
    - For macOS, install XCode
    - For Ubuntu, install a compiler through apt
2. GNU Scientific Library (GSL) version 1.0 or higher:
    - For macOS, install via Homebrew
    - For Ubuntu, install via apt
3. To compile using make:
    - Navigate to the main directory
    - Run `make serial`
    - The executable will be placed in the `./bin` directory

To build parallel NERDSS, you need to checkout to the `mpi` branch:

```
git checkout mpi
```

1. MPI:
    - For macOS, install OpenMPI with Homebrew: `brew install open-mpi`
    - For Ubuntu, install OpenMPI through apt: `sudo apt install openmpi-bin libopenmpi-dev`
2. To compile using make:
    - Navigate to the main directory
    - Run `make mpi` (with profiling support: *make mpi ENABLE_PROFILING=1*)
    - The executable will be placed in the ./bin directory

### Run Simulations

#### Run a quick trial with Google Colab

Click the following link to make a copy of the iPython notebook in your Google Colab and following the instructions on the Notebook (Note currently the link is under NERDSS development repo. The link will need to be updated once synced to released repo)

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mjohn218/nerdss_development/blob/master/docs/Run_NERDSS_colab.ipynb?copy=true)

#### Compile and run NERDSS on your local machine

1. Example input files are located in the subdirectories within the `sample_inputs` folder.

2. To start the serial simulation, use the command `./nerdss -f parms.inp`.

3. To start the parallel simulation, use the command `mpirun -np $nprocs  ./nerdss_mpi -f parms.inp`.

4. To debug the parallel code, use the command `mpirun -np 2 xterm -e gdb --ex 'b error' --ex r --args ./nerdss_mpi -f parms.inp -s 123`.

### Analyzing Results

1. Use the ioNERDSS PyPi library for visualizing simulation results.
2. Install ioNERDSS with *pip install ioNERDSS*.
3. Refer to the [ionerdss repository](https://github.com/JohnsonBiophysicsLab/ionerdss) for more details.

### Style Guides for Developers

1. We use the Google C++ Style Guide and the C++ Core Guidelines and prefer the C++ Standard Library and modern C++ features.
2. Provide comments and documentation to explain complex or non-obvious code sections. Write unit tests to ensure code quality.
