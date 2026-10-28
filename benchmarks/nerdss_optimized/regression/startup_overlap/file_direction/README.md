# The -c overlap pass, whichever molecule came from the file

With `-c`, `fixOverlappingMolecules()` compares each molecule the coordinate
file placed, always as the first molecule, with every other one, and re-places
the other one (a random molecule) on an overlap. `checkAndResolveOverlap()` used
to test only the first molecule's interfaces that react with the second's type,
each with the sigma of the *first* such reaction, against all of the second
molecule's interfaces. The full startup pass compares each pair both ways round
and so applied the union; the `-c` pass does not, so the same two molecules in
the same place overlapped or not depending on which one the file held. The
per-pair test (`interfacesOverlap()`, see `../per_pair`) is symmetric.

The model: one P (interface `p`) and one Q in a 100 nm box. Q has `q`, which
binds `p` at sigma 1 nm and is listed first; `s`, which binds `p` at sigma
2.5 nm; and `r`, which reacts with nothing.

`check.py` runs the model without `-c` to find where the two molecules start at
a seed (reading the file does not draw random numbers, so a `-c` run at the same
seed starts them in the same place), then with a one-molecule PDB file that
puts P, or Q, at a chosen distance from the other molecule's interface:

| placement        | overlap? | before: file P | before: file Q |
|------------------|----------|----------------|----------------|
| `p` 0.3 nm from `r` | no    | Q moved        | P stayed       |
| `p` 1.8 nm from `s` | yes   | Q stayed (`p` held `s` to `q`'s sigma 1) | P moved |
| `p` 0.3 nm from `q` | yes (control) | Q moved | P moved     |

Every other interface pair is outside its sigma in every case. The random
molecule must move exactly when the placement overlaps, both ways round.

```bash
./check.py <path-to>/nerdss
```

runs seeds 1-3 (`SEEDS` overrides), 21 runs of a few hundred ms. Before the fix
(on PR #91's c1f7ec32) it fails 6 checks, two per seed; with it, it passes.
