# Tests

Run everything with:

```
make test
```

from the repository root. This builds `ropebwt3` and `test/test_ps4g`, then runs both suites:

* **`test_ps4g.c`** — unit tests for `ps4g.c` (the gamete table, event
  accumulator, BED reader, and PS4G/`.npy` writers) in isolation, with no
  index building and no CLI. All expected values are hand-computed from a
  small synthetic set of events; see the comments in `test_finalize()`.

* **`run_integration.sh`** — builds the tiny pangenome index in
  `docs/examples/` (the same fixture used by `docs/examples/run.sh` and
  `docs/examples/refmap.out`) and runs `ropebwt3 refmap` end to end,
  including `--ps4g`, `--npy`, and `--label-bed`. Checks the PS4G file
  structure, the `#TotalUniqueCounts` invariant, that every `gameteSet` index
  is valid, the `.npy` header/shape (readable without a numpy dependency),
  and the diploid label columns. If a `python3` with `numpy` is importable
  (a few known conda envs are tried as a fallback) it also does a deeper
  `numpy.load` value check; if `valgrind` is installed it runs the tool
  under `--leak-check=full`. Both are skipped (not failed) when unavailable,
  with a `SKIP:` note on stderr.

Outputs land in `test/output/` and `test/tmp/` (gitignored); `make clean`
removes them along with the build artifacts.
