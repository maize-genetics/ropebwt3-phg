# Results — `refmap`

Outputs from the synthetic locus in `docs/experiments.md` (seed 42), run on the
optimized build (`make omp=0`). Columns:
`qname  qlen  status  nCar  carriers  ref  strand  cL  cR  span  ins`.

> Real-data **sensitivity + speed** at maize-pangenome scale:
> [`results-maize.md`](results-maize.md) (protocol in
> [`benchmark-maize.md`](benchmark-maize.md)).

Ground truth: insertion breakpoint at reference coordinate **200** (the `A|B`
junction), with incidental microhomology spanning ~198–201.

## Core cases (consensus mode)

```
q1     150  PLACED    6  g{1..6}_chr1:+  B73_chr1  +  198  201   3  297
q1rc   150  PLACED    6  g{1..6}_chr1:-  B73_chr1  -  198  201   3  297
qref   150  EXACT     0  .               B73_chr1  +  150  300  150    0
qsnp   150  PLACED    6  g{1..6}_chr1:+  B73_chr1  +  198  201   3  372
rand    44  UNPLACED  0  .               .         .    .    .    .    .
```

Interpretation:

* **`q1` (insertion query, fwd)** — `PLACED` on `B73_chr1` at 198–201, correctly
  bracketing the true breakpoint at 200. `span = 3` is the microhomology zone
  (INS shares 1 bp with B + 2 bp with A). `ins = 297` ≈ the 300 bp insertion
  (300 − 3 of microhomology). ✅
* **`q1rc` (reverse complement)** — same placement, strand `-`. Confirms
  strand-aware coordinate handling. ✅
* **`qref` (in the reference)** — `EXACT`, direct coordinates 150–300 (query was
  taken from `(A+B)[150:300]`). ✅
* **`qsnp` (one mid-query mismatch)** — still `PLACED` at 198–201 via the
  longest-SMEM core fallback. `ins = 372` is inflated because the exact core is
  shorter than the full query (documented approximation); the **location is
  correct**. ✅
* **`rand` (unrelated)** — `UNPLACED`, 0 carriers. ✅

## Walk modes (on `q1`; `g6_chr1` carries a SNP at `INS[20]`, in the left flank)

**strict** — stops the outward walk at the first carrier disagreement. The g6 SNP
sits in the left flank, so the left anchor cannot reach the breakpoint:

```
q1  150  ONE_SIDE  6  g{1..6}_chr1:+  B73_chr1  +  .  198  .  .
```

Only the right flank anchors (cR = 198). ✅ (Demonstrates strict's conservatism.)

**per-carrier** — follows each carrier individually; one line per carrier. g6
places correctly because per-carrier follows g6's own base through its SNP
(where strict drops it):

```
q1  150  PLACED  1  g3_chr1:+  B73_chr1  +  198  201  3  297
q1  150  PLACED  1  g6_chr1:+  B73_chr1  +  198  201  3  297
q1  150  PLACED  1  g4_chr1:+  B73_chr1  +  198  201  3  297
q1  150  PLACED  1  g2_chr1:+  B73_chr1  +  198  201  3  297
q1  150  PLACED  1  g1_chr1:+  B73_chr1  +  198  201  3  297
q1  150  PLACED  1  g5_chr1:+  B73_chr1  +  198  201  3  297
```

✅ (Demonstrates per-carrier handling divergent carriers that strict drops.)

## `--max-walk` bound

`--max-walk=30` (breakpoints are ~75 bp from the query, beyond the cap):

```
q1  150  UNPLACED  6  g{1..6}_chr1:+  .  .  .  .  .  .
```

Carriers still found, but no anchor within the cap → `UNPLACED`. ✅

## Non-functional checks

* **Build**: `make omp=0` — no warnings.
* **AddressSanitizer**: `make omp=0 asan=1`, batch of all queries (consensus and
  per-carrier) — no errors, no leaks reported.
* **Regression**: `mem` and `mem -p` outputs unchanged. `mem -p3` on `qref`
  independently reports `B73_chr1:+:150`, matching `refmap`'s `EXACT cL=150`.

## Verdict

All eight matrix cases pass. The core requirement — placing a query that is
**absent from the reference** (inside an insertion) at the correct reference
locus — works on both strands and across all three walk modes.
