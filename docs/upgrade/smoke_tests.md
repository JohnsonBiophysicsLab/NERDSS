# Serial Smoke Test Runner

Agent B owns the Phase 0 baseline build and smoke runner. The runner builds
serial NERDSS and executes a deliberately tiny generated input in an isolated
temporary directory, then records the command metadata needed for local and CI
triage.

## Command

```bash
python3 tools/run_smoke_tests.py --artifact-dir /tmp/nerdss-smoke
```

By default, the runner:

- checks that `gsl-config`, `make`, and a C++ compiler are available;
- runs `make serial`;
- stages `smoke.inp` and `A.mol` into a temporary directory;
- runs `bin/nerdss -f smoke.inp -s 123`;
- captures stdout, stderr, exit code, runtime, and produced file names/sizes;
- writes `report.json`, command logs, and the isolated run directory when
  `--artifact-dir` is provided.

The smoke input is intentionally minimal and is not a scientific regression
baseline. It only verifies that the serial executable can build, parse an input,
initialize one molecule, run a few iterations, and emit normal output files.

## Useful Options

```bash
python3 tools/run_smoke_tests.py --json
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss
python3 tools/run_smoke_tests.py --build-system cmake --artifact-dir /tmp/nerdss-smoke-cmake
```

If GSL or the selected build tool is missing, the runner exits with status `2`
and reports a skipped smoke with an actionable dependency message instead of
failing deep inside the build.
`--skip-build` bypasses build-tool dependency requirements and uses the provided
executable directly.

## Artifact Layout

When `--artifact-dir` is supplied, the directory contains:

- `report.json`: structured summary for CI or baseline archiving;
- `build.stdout.txt` and `build.stderr.txt`;
- `run.stdout.txt` and `run.stderr.txt`;
- `smoke_run/`: staged input files and all files produced by NERDSS.

Generated build outputs remain under the normal ignored build directories
(`bin/`, `obj/`, or `build/smoke-serial/`) and should not be committed.
