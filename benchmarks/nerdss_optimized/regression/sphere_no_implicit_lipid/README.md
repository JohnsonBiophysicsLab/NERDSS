# Spherical boundary without an implicit lipid

This model exists for the RS3D look-up crash. Every spherical model under
`sample_inputs/` carries an implicit lipid, so nothing in either case table
runs a sphere whose `RS3Dvect` table was never built, and the crash survived
in the sphere sweeps for as long as it did.

It is `sample_inputs/VALIDATE_SUITE/trimer` with two lines changed: the
`WaterBox = [118.41,118.41,118.41]` line becomes `isSphere = true` /
`sphereR = 73.46`, and `nItr` drops to 8000. The three `.mol` files are
byte-identical to the trimer's. 73.46 nm is the radius of equal volume --
1660510 nm^3 against the box's 1660218, 0.02% apart -- so the same model runs
at the same concentration in either geometry and the two can be compared.

Run it with:

```bash
mkdir -p /tmp/sphere_no_il && cp *.mol parms.inp /tmp/sphere_no_il && cd /tmp/sphere_no_il
<path-to>/nerdss -f parms.inp -s 20260810
```

## What it is built to show

Before the fix, the `-O3` build dies on the first step:

```
*************** BEGIN SIMULATION ****************
Segmentation fault                  (exit 139)
```

and the ASan build names the line:

```
ERROR: AddressSanitizer: SEGV on unknown address 0x000000000c80
The signal is caused by a READ memory access.
Hint: address points to the zero page.
    #0 sweep_separation_complex_rot_sphere.cpp:63
    #1 sweep_separation_dispatch.cpp:25
```

`0xc80` is 3200, which is `400 * sizeof(double)`: `RS3Dvect[RS3Dindex + 400]`
read through a null data pointer, because
`initialize_paramters_for_implicitlipid_and_compartment_model()` fills the
table only `if (systemIL == true || membraneObject.hasCompartment == true)`.

After the fix the same command exits 0 over 8000 iterations, and the ASan
build reports no error once `countTransition` is switched off in the three
`.mol` files -- with it on, the run stops instead on the unrelated
container-overflow in `associate_sphere.cpp:456` that reads a destroyed
complex's cleared `numEachMol`. That one is not specific to spheres: the
stock box trimer, unmodified, trips the identical overflow at the box twin
`associate_box.cpp:914`.

## Sanity of the result

Containment holds. Over seeds 11111, 20260810, 22222, 33333 and 44444, the
largest distance from the origin of any of the 900 printed points in
`DATA/final_coords.xyz` is 73.429, 73.443, 73.433, 73.410 and 73.277 nm,
every one inside `sphereR = 73.46`, and none of the 4500 points is outside.
The radial distribution is what a filled sphere gives: a pooled median of
56.75 nm against the 58.31 of a uniform fill, pulled slightly inward because
the boundary holds a whole complex inside, not just its centre.

The chemistry matches the equal-volume box. Summing the three bound-pair
columns of `DATA/copy_numbers_time.dat` at the last written step, over the
same five seeds in that order:

| geometry | per seed | mean | sd | sem |
| --- | --- | --- | --- | --- |
| box | 25, 29, 25, 31, 32 | 28.4 | 3.3 | 1.5 |
| sphere | 20, 20, 28, 29, 36 | 26.6 | 6.8 | 3.0 |

a difference of 1.8 +/- 3.4, consistent with zero. Five seeds is enough to
say the sphere is not obviously wrong and not enough to say the two agree to
better than about 10%; nothing here is a rate validation.
