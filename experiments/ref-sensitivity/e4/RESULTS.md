# E4 prototype — synteny-prior placement via carrier→B73 liftover

Goal: replace refmap's slow outward **walk** with a coordinate **projection**.
A carrier hit's exact coordinate comes free from the SSA at locate time; project
it to an approximate B73 coordinate through a precomputed carrier→B73 liftover,
and mark **NULL** (→ UNPLACED) wherever no confident projection exists, to keep
precision (per the "two SSA + null slots" design).

This is a **Python prototype** over `ropebwt3 mem -p` output (not yet the
integrated C "second SSA"), built only from the existing index — no external
aligner, no synteny data.

## Pipeline
1. `gen_anchors.py` — sample anchor windows along each carrier (stride S, len L).
2. `mem -p` the anchors → for every single-copy-per-taxon anchor with one B73
   hit, a pair (carrier_pos ↔ B73_pos).
3. `liftover.py` — sort pairs per (carrier,chr); `project()` interpolates between
   flanking anchors when they agree on chr and are locally collinear (slope≈1),
   else returns NULL.
4. `place_reads.py` — EXACT if a B73 hit; else project carrier hits, place at the
   agreed chr/median, else UNPLACED; MULTI if occ > max_occ.

## Liftover self-test (leave-one-out, stride 10k, L=100, max_occ=4)
- 179,517 usable liftover points from 639,501 anchors (28.1%).
- Projection: 85.7% committed, **14.3% NULL**; of committed, **99.76% within 5 Mb**.
- |Δ B73| median **7 bp**, p90 5.3 kb, p99 44 kb. Interpolation is very accurate.

## Read placement vs refmap (100k reads, tol ±5 Mb)

| variant | precision | recall | wrong-chr | placement speed |
|---|---:|---:|---:|---:|
| E0 baseline refmap | 57.3% | 54.6% | 35,768 | 40 s (walk) |
| E3 two-flank+max-occ=4 | 96.0% | 45.5% | 920 | 37 s (walk) |
| **E4 liftover (stride 10k, occ4)** | **95.9%** | **63.7%** | 1,068 | ~3 s locate + project |

By status: EXACT 42,242 (96.7% correct) · PLACED 24,218 (**94.5% correct**) ·
MULTI 18,062 (filtered) · UNPLACED 15,045 (null/disagree) · AMBIG 433.

**Takeaway:** the synteny prior beats the walk on precision, recall, AND speed at
once. The carrier-projected PLACED bucket is 94.5% correct vs refmap's 14.7%
ONE_SIDE. Recall +18 pts over E3 at the same ~96% precision, ~15× faster.

## Parameter experiments
(running — denser anchors, E4 max_occ sweep, agree-tol)
