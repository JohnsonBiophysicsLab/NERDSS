# Coordinate file formats for `-c`

`-c <file>` used to read PDB only. `src/io/read_coordinate_file.cpp` now reads
PDB, XYZ and CSV, picks the format from the extension or, failing that, from
the content, and `updateMoleculeCoordinates()` in
`src/system_setup/generate_coordinates.cpp` consumes the same name -> COM map
it built itself before. The formats, the detection rules and the coordinate
frame are documented in `docs/mainpage.md` (section "Initial coordinates from
a file") and in `include/io/read_coordinate_file.hpp`.

`check.py` runs a three-type model (P, Q, Longname: one interface each, 0.5 nm
from the COM, reaction P(p)+Q(q), `nItr = 1`) and reads back
`DATA/initial_crds.xyz`, which a one-step run leaves as the whole result:

1. **Equivalence.** The same placements written as PDB, XYZ, CSV with a
   header and CSV without one, run with the same seed, leave byte-identical
   `initial_crds.xyz`; the placed molecules sit at the file positions and the
   others keep the positions of a run without `-c` (the RNG stream is the
   same up to the moment the file is read).
2. Names longer than PDB's 4 characters work in XYZ and CSV.
3. **Detection.** Upper/mixed-case extensions; no or an unknown extension for
   each format (content sniffing, including a `HETATM`-only PDB and a
   headerless `.dat` CSV); a file no rule matches stops the run.
4. **Edge cases.** XYZ with a text or empty comment line, CRLF endings in XYZ
   and CSV, spaces around CSV fields and blank lines, trailing blank lines, a
   CSV whose first row is numeric (not a header), a header-only CSV, an
   unknown molecule name (skipped), an XYZ with more rows than line 1
   announces (warning, extras ignored) or fewer (error), `0` rows.
5. **Malformed rows** stop the run with the file name and the 1-based line
   number: a non-number, a numeric prefix (`5.5x`), a wrong field count, an
   empty name, a line 1 that is not a whole number. A missing file keeps the
   old behaviour (message, generated coordinates kept, exit 0).
6. With a `nerdss_mpi` path: the four formats at np=1 and np=2 agree with each
   other and a malformed row stops the run there too.

```
python3 benchmarks/nerdss_optimized/regression/coordinate_file_formats/check.py bin/nerdss [bin/nerdss_mpi] [work-dir]
```

Runs land in the work dir (a fresh temp dir by default); the script prints one
PASS/FAIL line per check and exits non-zero on any FAIL.
