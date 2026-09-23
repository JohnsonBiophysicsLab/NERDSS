# Restart format round trip

Restart files are JSON (the port of nerdss_development's `json-restarts`
branch); the positional `.dat` format that came before is still read, and
`LEGACY_write_restart()` can still write it.  `check.sh` tests that the two
formats carry the same state and that each reads back everything it writes:

```bash
make serial
./check.sh ../../../../bin/nerdss [work-dir]
```

`restart_convert.cpp` reads a restart file of either format (`read_restart()`
tells them apart by the first byte) and writes it back in both.  It is built
against `obj/serial`, so the reader and writers under test are the binary's
own.  For each case in `cases.tsv` the model is run, and then

```
J1 = DATA/restart.json, written at the end of the run
restart_convert J1  ->  J2 (JSON) and D1 (.dat)
restart_convert D1  ->  J3 (JSON) and D2 (.dat)
```

must satisfy

| comparison | shows |
| --- | --- |
| `J1 == J2`, byte for byte | the JSON reader consumes everything the JSON writer wrote |
| `D1 == D2`, byte for byte | the `.dat` reader consumes everything the `.dat` writer wrote |
| `J1 == J3` on every value but `stateChangeIface`, within 1e-20 | a `.dat` file holds the state a JSON file holds; `stateChangeIface` is the one record the `.dat` format never carried |

The 1e-20 is what the `.dat` format keeps: it prints most doubles in fixed
notation with twenty decimal places, so a template interface coordinate of
4.7e-15 comes back from it with six significant digits, and a lifetime of
2.41e-05 loses its last bit.  JSON prints every double with the digits it needs
to read back exactly, which is why `J1 == J2` is asked for byte for byte.

With `REFERENCE_BINARY` set to a build that still writes `DATA/restart.dat`,
the same model is also run with it and its final `restart.dat` must equal `D1`
byte for byte: the JSON file, converted, reproduces what the old build wrote.

The cases are chosen so that every section of the file is exercised at least
once; the table says what each one adds.  Exit status 0 means every case
passed.
