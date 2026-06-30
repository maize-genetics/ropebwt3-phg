# `ropebwt3 refmap` — usage guide

Place a query sequence on a designated **reference genome** inside a multi-genome
(pangenome) index — including when the query is **absent from the reference**
because it falls inside an insertion (a deletion on the reference side) and so
has no exact match there.

A complete, runnable example lives in [`docs/examples/`](examples/); this guide
walks through it.

---

## 1. When to use it

You have one BWT index built from many genomes of a species (e.g. the 26 maize
NAM founders), and one of them is the reference (e.g. B73). You have a short
query (say 150 bp) that hits several genomes but **not** the reference, because
the region is part of an insertion relative to the reference. `mem` and `sw`
cannot give you a reference coordinate for it; `refmap` can, by walking outward
through the carrier genomes to the insertion breakpoints and re-anchoring there.

If the query *is* present in the reference, `refmap` reports its exact
coordinates (`status = EXACT`) instead.

---

## 2. Requirements

`refmap` needs the full index triad next to `<base>`:

| File | Built by | Purpose |
|---|---|---|
| `<base>.fmd` | `ropebwt3 build -d` | the BWT itself |
| `<base>.fmd.ssa` | `ropebwt3 ssa` | sampled suffix array → coordinates |
| `<base>.fmd.len.gz` | `seqtk comp` / awk | sequence names and lengths |

The reference is selected by **sequence-name prefix**, so name the reference
sequences with a common prefix (e.g. `B73_chr1`, `B73_chr2`, …) distinct from
the other genomes.

On macOS build the tool with `make omp=0` (Apple clang lacks `-fopenmp`).

---

## 3. Input files

### Pangenome FASTA — [`docs/examples/pangenome.fa`](examples/pangenome.fa)

Seven sequences: the reference `B73_chr1` (= flank `A` + flank `B`) and six
carriers (= `A` + a 500 bp insertion + `B`). One carrier (`Ki3_chr1`) has a SNP
inside the insertion.

```
>B73_chr1
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG...
>Mo17_chr1
...
>Ki3_chr1
...
>B97_chr1
...
```

### Query FASTA — [`docs/examples/queries.fa`](examples/queries.fa)

```
>ins_query 150 bp from the middle of the insertion; absent from the reference
GTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAGCCGGTGACTCCTAATGCTAA...
>ins_query_rc reverse complement of ins_query
...
>ref_query spans the A|B junction; present in the reference
...
>ins_query_1snp ins_query with one mismatch (no end-to-end carrier match)
...
>unrelated random sequence, in no genome
...
```

---

## 4. Building the index and running

The script [`docs/examples/run.sh`](examples/run.sh) does all of it. From the
repository root, after `make omp=0`:

```sh
sh docs/examples/run.sh
```

which runs:

```sh
# 1. BWT (static FMD format)
ropebwt3 build -d -o docs/examples/pangenome.fmd docs/examples/pangenome.fa

# 2. sampled suffix array
ropebwt3 ssa -s8 -o docs/examples/pangenome.fmd.ssa docs/examples/pangenome.fmd

# 3. names + lengths (seqtk alternative shown; the script uses awk)
seqtk comp docs/examples/pangenome.fa | cut -f1,2 | gzip > docs/examples/pangenome.fmd.len.gz

# 4. place the queries on the B73 reference
ropebwt3 refmap --ref-prefix=B73 docs/examples/pangenome.fmd docs/examples/queries.fa
```

---

## 5. Output

One tab-separated line per query. Full captured output (with a column legend and
the strict / per-carrier variants) is in
[`docs/examples/refmap.out`](examples/refmap.out).

```
qname          qlen status   nCar carriers                       refName  strand cL  cR  span ins
ins_query       150 PLACED      6 Oh43,B97,CML247,Mo17,Ki3,Tx303 B73_chr1 +      300 300 0    500
ins_query_rc    150 PLACED      6 Oh43,B97,CML247,Mo17,Ki3,Tx303 B73_chr1 -      300 300 0    500
ref_query       150 EXACT       0 .                              B73_chr1 +      250 400 150  0
ins_query_1snp  150 PLACED      6 Oh43,B97,CML247,Mo17,Ki3,Tx303 B73_chr1 +      300 300 0    575
unrelated       120 UNPLACED    0 .                              .        .      .   .   .    .
```

### Columns

| # | Field | Meaning |
|---|---|---|
| 1 | qname | query name |
| 2 | qlen | query length |
| 3 | status | `PLACED`, `ONE_SIDE`, `EXACT`, or `UNPLACED` |
| 4 | nCarrier | number of carrier genomes the query was found in |
| 5 | carriers | carrier names with strand (`name:+`/`name:-`), or `.` |
| 6 | refName | reference sequence placed on, or `.` |
| 7 | strand | reference strand (`+`/`-`) |
| 8 | cL | reference coordinate of the **left** breakpoint, or `.` |
| 9 | cR | reference coordinate of the **right** breakpoint, or `.` |
| 10 | refSpan | `cR - cL`; reference bases spanned (0 = clean insertion point) |
| 11 | insSize | implied inserted size in the carriers (bp) |

### How to read the example

* **`ins_query`** is absent from B73 but found in all six carriers. `refmap`
  places it at **B73_chr1:300**, exactly the `A|B` junction, with `refSpan = 0`
  (a clean insertion point) and `insSize = 500` (the true insertion length). This
  is the headline use case: a query that *cannot* be found in the reference by
  exact matching gets a correct reference coordinate.
* **`ins_query_rc`** gives the same locus on the `-` strand — coordinates are
  strand-aware.
* **`ref_query`** is in the reference, so `status = EXACT` and `cL`/`cR` are its
  direct coordinates (250–400).
* **`ins_query_1snp`** has a mismatch and matches no carrier end to end, yet is
  still placed at 300 via its longest exact core. Note `insSize = 575` is
  **approximate** in this fallback case (the exact core is shorter than the
  query); the location is correct.
* **`unrelated`** is in no genome → `UNPLACED`.

### Status values

| status | meaning |
|---|---|
| `PLACED` | both flanks re-anchored; `cL`/`cR` bracket the query on the reference |
| `ONE_SIDE` | only one flank anchored (the other side diverges or is too far) |
| `EXACT` | the query occurs in the reference itself; `cL`/`cR` are its coordinates |
| `UNPLACED` | the query is in no genome, or no anchor within `--max-walk` |

> **Microhomology.** With real sequence, the inserted bases often share a few
> bases with the flanks, so `cL` and `cR` can differ by a few bp and `refSpan`
> is small but non-zero. That is correct — it reflects genuine breakpoint
> ambiguity, not an error. (This example uses sequence with no microhomology, so
> `cL == cR == 300`.)

---

## 6. Options

```
--ref-prefix=STR   reference = sequences whose name starts with STR   [required]
--max-walk=NUM     max bases to walk outward along carriers per flank [5000]
--walk-mode=STR    carrier path: consensus | strict | per-carrier     [consensus]
-l INT             min anchor length when re-mapping a flank          [19]
-t INT             number of threads                                  [4]
-L                 one sequence per line in the input
```

### `--walk-mode`

When carriers disagree within a flank (here `Ki3_chr1` has a SNP in the
insertion):

* **`consensus`** (default) follows the base supported by the most carriers, so a
  single divergent carrier does not derail the walk.
* **`strict`** stops at the first disagreement. In the example the Ki3 SNP sits in
  the left flank, so the left anchor fails and `ins_query` becomes `ONE_SIDE`
  (only `cR` reported).
* **`per-carrier`** walks each carrier individually and prints one line per
  carrier — useful when carriers carry different insertion alleles. Ki3 still
  places correctly because its own path is followed through its SNP.

See the strict and per-carrier blocks in
[`docs/examples/refmap.out`](examples/refmap.out).

### `--max-walk`

Caps how far each flank walks outward. A query buried inside an insertion larger
than the cap stays `UNPLACED` (no reference anchor is reached). Raise it for
large structural variants at the cost of more work per query.
