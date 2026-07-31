# Tests

Run everything with:

```
make test
```

from the repository root. This builds `ropebwt3`, `test/test_ps4g`, and
`test/test_hitcount`, then runs all three:

* **`test_ps4g.c`** — unit tests for `ps4g.c` (the gamete table, event
  accumulator, BED reader, and PS4G/`.npy` writers) in isolation, with no
  index building and no CLI. All expected values are hand-computed from a
  small synthetic set of events; see the comments in `test_finalize()`. This
  includes the case that matters most: two distinct gameteSets sharing one
  bin must produce two separate `.npy` rows, never summed into one, and
  `--npy-binary` must clip counts to presence without touching the labels.

* **`test_hitcount.c`** — unit tests for `hitcount.c`, the atomic counter
  behind `--target-hits`, in isolation. Covers `init`/`add`/`reached`/`short`
  semantics, that `target<=0` means "always off" regardless of how many hits
  are added, and a real `pthread_create` concurrency test (several threads
  incrementing at once) proving the atomic add never loses an update — the
  actual failure mode this module exists to prevent, since the read stage and
  write stage run on different threads in production.

* **`run_integration.sh`** — builds the tiny pangenome index in
  `docs/examples/` (the same fixture used by `docs/examples/run.sh` and
  `docs/examples/refmap.out`) and runs `ropebwt3 refmap` end to end,
  including `--ps4g`, `--npy`, `--npy-binary`, `--label-bed`, and
  `--target-hits`. Checks the PS4G file structure, the `#TotalUniqueCounts`
  invariant, that every `gameteSet` index is valid, the `.npy` header/shape
  (readable without a numpy dependency), a count-vs-binary cell comparison,
  the diploid label columns, and (for `--target-hits`) an exact,
  hand-computed record count on a purpose-built 30-record fixture alternating
  known-EXACT and known-UNPLACED reads, plus a shortfall/warning case. If a
  `python3` with `numpy` is importable (a few known conda envs are tried as a
  fallback) it also does a deeper `numpy.load` value check; if `valgrind` is
  installed it runs the tool under `--leak-check=full` (including once with
  `--target-hits` under real multi-threaded compute). Both are skipped (not
  failed) when unavailable, with a `SKIP:` note on stderr.

Outputs land in `test/output/` and `test/tmp/` (gitignored); `make clean`
removes them along with the build artifacts.
