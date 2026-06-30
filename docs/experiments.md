# Experiments — `refmap`

No multi-genome pangenome data is available locally, so testing uses a small
**synthetic insertion locus** where the ground truth is known exactly. Results
are in `docs/results.md`.

## Synthetic design

Seven sequences sharing two flanks `A` (200 bp) and `B` (200 bp):

* `B73_chr1` = `A + B` — the **reference** (the insertion is deleted here).
* `g1_chr1 … g6_chr1` = `A + INS + B`, with `INS` a 300 bp insertion. The
  carriers share `INS`, except `g6_chr1` carries one SNP inside `INS` (at
  `INS[20]`) to exercise `consensus` vs `strict` vs `per-carrier`.

Queries:

* `q1` = `INS[75:225]` — 150 bp from the **middle of the insertion**; absent from
  the reference. Ground-truth breakpoint at reference coordinate **200**
  (the `A|B` junction).
* `q1rc` = reverse complement of `q1` — same locus, `-` strand.
* `qref` = `(A+B)[150:300]` — a query that **is** in the reference (positions
  150–300), to test the `EXACT` path.
* `qsnp` = `q1` with one mid-query mismatch — no end-to-end carrier match; tests
  the longest-SMEM fallback.
* `rand` = unrelated random sequence — tests `UNPLACED` with 0 carriers.

The test data has incidental **microhomology**: `INS` shares 1 leading base with
`B` and 2 trailing bases with `A`, so the exact breakpoint is ambiguous over
reference coordinates ~198–201. This is realistic and the expected output
reflects it (`cL`/`cR` bracket the zone).

## Reproduce

### 1. Build the tool

```sh
make omp=0
```

### 2. Generate the synthetic locus

```python
# gen_test.py
import random, gzip
random.seed(42)
rnd = lambda n: ''.join(random.choice('ACGT') for _ in range(n))
A, B, INS = rnd(200), rnd(200), rnd(300)
ref, car = A + B, A + INS + B
ins6 = list(INS); ins6[20] = 'A' if INS[20] != 'A' else 'C'
car6 = A + ''.join(ins6) + B
seqs = [("B73_chr1", ref)] + [(f"g{i}_chr1", car) for i in range(1, 6)] + [("g6_chr1", car6)]
open("test.fa", "w").write(''.join(f">{n}\n{s}\n" for n, s in seqs))
with gzip.open("test.fa.len.gz", "wt") as f:
    f.write(''.join(f"{n}\t{len(s)}\n" for n, s in seqs))
q = INS[75:225]
open("q.fa", "w").write(f">q1\n{q}\n")
comp = {'A':'T','C':'G','G':'C','T':'A'}
open("qrc.fa", "w").write(">q1rc\n" + ''.join(comp[c] for c in reversed(q)) + "\n")
open("qref.fa", "w").write(">qref\n" + (A+B)[150:300] + "\n")
qs = list(q); qs[75] = 'A' if qs[75] != 'A' else 'C'
open("qsnp.fa", "w").write(">qsnp\n" + ''.join(qs) + "\n")
open("rand.fa", "w").write(">rand\nGCGCTTAACCGGTTAACCGGAATTCCGGAATTCCAAGGTTCCAA\n")
```

```sh
python3 gen_test.py
```

### 3. Build the index (the full triad: .fmd + .ssa + .len.gz)

```sh
ropebwt3 build -d -o t.fmd test.fa
ropebwt3 ssa  -s4 -o t.fmd.ssa t.fmd
cp test.fa.len.gz t.fmd.len.gz
```

### 4. Run

```sh
ropebwt3 refmap --ref-prefix=B73 t.fmd q.fa      # insertion query (consensus)
ropebwt3 refmap --ref-prefix=B73 t.fmd qrc.fa    # reverse strand
ropebwt3 refmap --ref-prefix=B73 t.fmd qref.fa   # in-reference -> EXACT
ropebwt3 refmap --ref-prefix=B73 t.fmd qsnp.fa   # SMEM-core fallback
ropebwt3 refmap --ref-prefix=B73 t.fmd rand.fa   # UNPLACED
ropebwt3 refmap --ref-prefix=B73 --walk-mode=strict      t.fmd q.fa
ropebwt3 refmap --ref-prefix=B73 --walk-mode=per-carrier t.fmd q.fa
ropebwt3 refmap --ref-prefix=B73 --max-walk=30           t.fmd q.fa  # too short -> UNPLACED
```

## Test matrix

| Case | Input | Mode | Expectation |
|---|---|---|---|
| Insertion, fwd | `q1` | consensus | `PLACED` B73_chr1 `+` ~200 |
| Insertion, rev | `q1rc` | consensus | `PLACED` B73_chr1 `-` ~200 |
| In reference | `qref` | consensus | `EXACT` 150–300 |
| Mismatched query | `qsnp` | consensus | `PLACED` (longest core); insSize approx |
| Unrelated | `rand` | consensus | `UNPLACED`, 0 carriers |
| Divergent carrier | `q1` | strict | `ONE_SIDE` (left flank stops at g6 SNP) |
| Divergent carrier | `q1` | per-carrier | one `PLACED` line per carrier incl. g6 |
| Walk too short | `q1` | `--max-walk=30` | `UNPLACED` |

## Validation performed

* Builds clean with `make omp=0` (no warnings).
* Runs clean under `make omp=0 asan=1` (AddressSanitizer, batch of all queries).
* `mem` / `mem -p` regression-checked (unchanged; `mem -p` independently
  confirms the reference coordinate that `refmap` reports as `EXACT`).
