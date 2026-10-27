# Re-placed molecules after `-c` are checked again

With `-c coords.pdb`, `generate_coordinates()` first places every molecule at
random and resolves overlaps among them, then moves the first N molecules of
each type onto the file's positions and runs
`fixOverlappingMolecules(..., checkOnlyChanged = true)`. That pass compares each
file molecule with every other molecule; a random molecule that collides with a
file molecule is re-placed at random (`create_random_coords`), and two file
molecules that collide only print a WARNING.

Before the fix the re-placed molecule was never looked at again. It could land
on another random molecule, and the run started with that overlap while the
log's last pass said `|-Overlapping Proteins (false): 0`. Now every re-placed
molecule is also compared with all other molecules, and re-placed again while
it overlaps one, until the pass is clean or reaches its iteration cap. Only the
re-placed molecule ever moves: file molecules stay where the file put them, and
random molecules that collided with nothing keep their positions.

## Model

`P` and `Q`, one interface each 0.5 nm from the center, `P(p) + Q(q)` at
`sigma = 2`, `nItr = 1`. `parms.template` is filled in by `check.py`.

A run without `-c` and a run with `-c` use the same random stream up to the
moment the coordinate file is read, so a reference run gives the position of
every random molecule. In a 30 nm box with 120 P and 20 Q, the file puts the
first 20 P exactly on the 20 Q of the reference, so every Q has to be re-placed
among 100 random P.

For every seed `check.py` requires: no P/Q pair closer than sigma; the 20 file
P exactly at their file positions; the 100 random P exactly at their reference
positions; all 20 Q moved; and the last `|-Overlapping Proteins` line `false`.
A second case puts a file P and a file Q on the same spot: it must stay a
WARNING, not counted as an overlap, with neither molecule moved.

```bash
python3 -I check.py <path-to>/nerdss [work-dir]
SEEDS=300 python3 -I check.py <path-to>/nerdss
```

For `nerdss_mpi`, pass a wrapper that runs `mpirun -np N nerdss_mpi "$@"`;
each rank's log is read from `DATA/output_<rank>.txt`.

## Results

| build | seeds | seeds with overlaps left | overlaps left | max per seed |
|---|---|---|---|---|
| before (4017ece8) | 30 | 28 | 77 | 6 |
| after | 30 | 0 | 0 | 0 |
| before (4017ece8) | 300 | 274 | 720 | 8 |
| after | 300 | 0 | 0 | 0 |

In both builds the file P and the random P never moved; only the Q did.
