# Sample-number-relative thresholds (audit)

Several caps were hardcoded for a few-founder simulator and silently truncate on a
real pangenome. They must scale with **N = number of samples/taxa** (distinct name
prefix before `_`). **Key subtlety:** the FMD holds *both strands*, so a sequence
present once per taxon at one locus has an FM **interval size ≈ 2N** (and yields ~2N
located occurrences). Interval/locate caps therefore need ~2N; carrier/gamete-set
caps need ~N. Ed's rule: scale with N but **ceiling at 256** (keep counts byte-sized).

Compute N once from the taxa count (search.c, in the is_ref block) and derive the rest.

| # | constant | file | caps | was | now | status |
|--|--|--|--|--|--|--|
| 1 | `chain_max_occ` | search.c | SMEM FM-interval (chain skips SMEMs above) | 5 | `min(2N,256)` auto | **DONE** |
| 2 | `max_pos` (chain) | search.c | # located occ/SMEM → gamete set | 64 | `>= chain_max_occ` | **DONE** |
| 3 | `RB3_RM_MAX_CARRIER` | search.c | carriers kept per placed read (refmap) | 64 | `min(N,256)`, dynamic arrays | TODO |
| 4 | `tmp[256]` gamete cap | ps4g.c | gametes per PS4G event | 256 | keep 256 (byte ceiling) | keep |
| 5 | k-mer locate `pos[64]`/cap | lift.c | carrier hits per k-mer (lift build) | 64 | `min(2N,256)` | TODO |

Already correct (the pattern): refmap `--max-occ` `<0 = auto = #taxa`.

## Impact by graph
- **NAM (N=27):** #1 was catastrophic — B97 leaf RNAseq placement 1.9%→40% (3') /
  1.2%→61% (full-length), source-recovery 67%→96% / 83%→99% after the fix. #2/#5
  borderline (2·27=54 < 64). Now fixed for chain.
- **282-panel WGS / PanAnd across-genera (N>64, some >256):** trips #2,#3,#5 (>64)
  and #4 (>256). #3/#5 (refmap + lift) still TODO — needed before those rows, or
  the true founder is silently truncated out of the PS4G (founder-dropout).

## Remaining work (#3, #5)
Replace fixed-size stack arrays (`rsids[64]`, `rposs[64]`, `v[64]`, `car_seen[64]`,
`carriers[64]` in refmap_place*/refmap_query; `pos[64]` in lift.c) with N-sized
(min(N or 2N, 256)) allocations, and raise the corresponding `1<<16`/`64` locate
caps consistently. Lower risk to do together with a small helper that returns the
byte-capped N-derived bound.
