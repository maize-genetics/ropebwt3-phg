# Tests

Run everything with:

```
make test
```

from the repository root. This builds `ropebwt3`, `test/test_ps4g`, `test/test_hitcount`,
`test/test_lift` and `test/test_extract`, then runs them all:

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

* **`test_extract.c`** — unit tests for coordinate-based extraction: builds an
  FM-index in memory from known random sequences (N runs, a 1-bp sequence, a
  repeat), then checks that every `isa` entry point LF-walks exactly `pos` steps
  to its own string start, extraction at/around every entry, sequence starts and
  ends, length 1 and whole sequences on both strands for several bucket sizes,
  the error codes, `.ssa` streaming vs in-memory, `.isa` round trip/validation,
  and hand-computed `rb3_lift_project_inv` (reference->assembly) projections.

* **`test_lift.c`** — unit tests for the lift reference-position index and the
  `--anchor-dist-npy` ternary/distance block (`rb3_lift_nearest_ref`,
  `rb3_lift_ternary_state`, `rb3_ps4g_npy_finalize` with a lift index), including
  the per-segment judgement of a spliced (`chain`) read's exon rows.

`run_integration.sh` section 6 runs `chain --npy --anchor-dist-npy` end to end on
the fixture plus a synthetic spliced read: the missing-`--lift` error, unchanged
stdout/PS4G/bins.tsv versus plain `--npy`, the `3*n_gametes+2` shape and layout
sidecar (identical to `refmap`'s), and a numpy check of the ternary definition.

`run_integration.sh` section 7 covers `isa`/`extract`: regions and TSV queries against
the FASTA, error reporting/exit status, `--clamp`, `-t4` vs `-t1`, `-w`, in-memory vs
file entry points, `--lift`/`--project-only` (incl. inverted and split founders), and
a scaled-down run of the shared simulated pangenome (`simulate_pangenome.py`, checked
by `check_extract_fixture.py`; skipped without numpy).

Outputs land in `test/output/` and `test/tmp/` (gitignored); `make clean`
removes them along with the build artifacts.
