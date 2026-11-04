| case | np | runs/failed | own | struct | ASan | conflicts | mass | max\|t\| 1st | max\|t\| 2nd (p) | column (2nd half) |
|---|---|---|---|---|---|---|---|---|---|---|
| clathrin | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | 2.23 | 1.69 (0.15) | clat(cd3!1).clat(cd3!1) mpi=1.2 ser=2.2 |
| clathrin | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | 0.61 | 1.08 (0.31) | clat(cd1) mpi=90.8 ser=89.8 |
| clathrin | 4 | 4/0 | 0 | 0 | 0 | 0 | 4bad | 86.57 | 164.77 (4.9e-07) | clat(td1) mpi=2.5 ser=100.0 |
| closed_homoTrimer | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 0.22 (0.83) | A(c) mpi=923.0 ser=921.4 |
| closed_homoTrimer | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 1.85 (0.098) | A(c) mpi=935.5 ser=921.4 |
| closed_homoTrimer | 4 | 4/0 | 0 | 0 | 0 | 0 | 4bad | - | 162.48 (1.3e-14) | A(c) mpi=5.0 ser=921.4 |
| hetTrimer | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 1.10 (0.34) | A(a1) mpi=222.8 ser=227.4 |
| hetTrimer | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 1.29 (0.24) | B(b1) mpi=232.5 ser=225.2 |
| hetTrimer | 4 | 4/0 | 0 | 0 | 0 | 0 | 4bad | - | 179.63 (2.7e-15) | B(b2) mpi=1.5 ser=227.4 |
| hexamer | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 0.27 (0.8) | hex(cd1) mpi=942.8 ser=943.9 |
| hexamer | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 0.10 (0.92) | hex(cd1) mpi=944.2 ser=943.9 |
| hexamer | 4 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 0.48 (0.64) | hex(cd1) mpi=946.0 ser=943.9 |
| homoTrimer | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 0.22 (0.83) | A(c) mpi=923.0 ser=921.4 |
| homoTrimer | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 1.85 (0.098) | A(c) mpi=935.5 ser=921.4 |
| homoTrimer | 4 | 4/0 | 0 | 0 | 0 | 0 | 4bad | - | 162.48 (1.3e-14) | A(c) mpi=5.0 ser=921.4 |
| mem_localization | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 3.02 (0.024) | B(m) mpi=98.6 ser=99.5 |
| mem_localization | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 2.50 (0.036) | M(m) mpi=3752.5 ser=3753.8 |
| mem_localization | 4 | 4/0 | 0 | 0 | 0 | 0 | ok | - | 0.33 (0.75) | M(m) mpi=3753.6 ser=3753.8 |
| michaelis_menten | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | 1.34 | 0.63 (0.55) | S(ser~P) mpi=0.7 ser=0.5 |
| michaelis_menten | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | 471.17 | 157.37 (1.1e-13) | S(ser~U) mpi=0.0 ser=103.7 |
| michaelis_menten | 4 | 4/0 | 0 | 0 | 0 | 0 | 4bad | 37.72 | 32.07 (2.6e-09) | S(ser~U) mpi=135.8 ser=103.7 |
| rev_2D | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | 0.79 | 0.64 (0.56) | A(a!1).R(r!1) mpi=114.8 ser=118.2 |
| rev_2D | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | 1.71 | 0.79 (0.47) | A(a) mpi=685.7 ser=681.8 |
| rev_2D | 4 | 4/0 | 0 | 0 | 0 | 0 | ok | 0.63 | 0.29 (0.79) | A(a!1).R(r!1) mpi=120.2 ser=118.2 |
| rev_3D | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | 0.07 | 0.74 (0.49) | A(a) mpi=549.2 ser=553.9 |
| rev_3D | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | 1.22 | 0.39 (0.71) | A(a!1).R(r!1) mpi=448.7 ser=446.1 |
| rev_3D | 4 | 4/0 | 0 | 0 | 0 | 0 | ok | 0.35 | 2.36 (0.042) | A(a) mpi=543.4 ser=553.9 |
| rev_3Dto2D | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | 1.06 | 1.58 (0.16) | M(r) mpi=2899.8 ser=2909.8 |
| rev_3Dto2D | 2 | 4/0 | 0 | 0 | 0 | 0 | ok | 3.26 | 2.74 (0.024) | A(a!1).M(r!1) mpi=104.2 ser=90.2 |
| rev_3Dto2D | 4 | 4/0 | 0 | 0 | 0 | 0 | ok | 0.57 | 1.07 (0.32) | M(r) mpi=2903.9 ser=2909.8 |
| trimer | 1 | 4/0 | 0 | 0 | 0 | 0 | ok | 2.14 | 2.71 (0.023) | C(c2!1).B(b1!1) mpi=9.3 ser=8.3 |
| trimer | 2 | 4/0 | 0 | 0 | 0 | 0 | 4bad | 53.78 | 46.96 (1.9e-11) | A(a2) mpi=50.2 ser=92.0 |
| trimer | 4 | 4/0 | 0 | 0 | 0 | 0 | 4bad | 32.13 | 26.08 (3.5e-09) | A(a2) mpi=68.5 ser=92.0 |
| unimol_state_change | 1 | 4/3 | 0 | 0 | 0 | 0 | ok | - | - |  |
| unimol_state_change | 2 | 4/3 | 0 | 0 | 0 | 0 | 1bad | - | - |  |
| unimol_state_change | 4 | 4/4 | 0 | 0 | 0 | 0 | ok | - | - |  |
