# Paired-end reads in ropebwt3 — plan and phase-1 status

![workflow](plan_diagram.png)

*Regenerate the diagram (`python3 docs/paired-end/plan_diagram.py`) whenever the stage list changes.*

## Scope and decisions (2026-10-14)

| Question | Decision |
|---|---|
| Base | `local/paired-end` = `local/placed-fix` + merge of `chain-founder-support` (one binary for refmap + chain) |
| Input | `--paired` + two files (R1 R2) **or** one interleaved file; mate names must match up to `/1` `/2` |
| Insert size | bwa-style estimate per batch from EXACT/EXACT pairs (≥25), or fixed `--pe-isize=MEAN,SD` |
| Output | one PS4G/npy event per pair; founder set = intersection of the mates' sets |
| Phasing | phase 1 = I/O + insert size + refmap; phase 2 = chain; phase 3 = mem/sw |
| Simulator | stays single-end for now (separate step) |

Without `--paired` nothing changes: single-end output is byte-identical to the parent binaries.

## Phase 1: refmap (implemented)

**Pair classes** (`rb3_pe_resolve`; printed as a trailing stdout column, paired mode only):

| Class | When | PS4G/npy |
|---|---|---|
| `CONC` | both mates placed on one contig, `|Δpos| ≤ max_ins + --pe-slack` | 1 event at the anchor (EXACT beats PLACED; tie → mate 1) |
| `RESC` / `MOVED` | discordant, but one mate has an alternative position near the other: a PLACED mate's other lift projection, or an EXACT mate's other B73 copy. The less confident mate is moved first | the `MOVED` mate gets that position; 1 event at the anchor |
| `SOLO` | exactly one mate placed | 1 event at it |
| `DISC` | discordant and not rescuable | each mate is its own single-end event |
| `NONE` | neither mate placed | none |

**Founder set of a pair event** (`rb3_pe_event_gametes`): the sorted intersection of the two
mates' sets. If the intersection is empty, or the mate has no set, the anchor's set is used.
If one mate's occurrence sample was capped (`RB3_RM_MAX_ASSEMBLY`), it cannot veto founders and
the other mate's set is used. A `MOVED` mate first keeps only the founders whose own
projection landed near the anchor (`rb3_pe_near_gametes`). Without that, a TE/paralog copy
in another founder would count as sharing the fragment.

**`--pe-trust-exact` (opt-in):** a discordant EXACT + PLACED pair with no alternative near
the other mate moves the PLACED mate onto the EXACT one. It is off by default because on
non-collinear (insertion/TE) fragments the EXACT mate is often at a repeat copy (results below).

**Defaults:**
- `--pe-slack=2000`. That is the lift-projection error scale, the same as `--anchor-dist-thresh`.
- `max_ins = max(mean + 4·sd, p75 + 3·IQR)`.
- 1000 bp is used until a batch has 25 EXACT/EXACT pairs, with a one-time warning.

**Not supported with `--paired` (exits with an error):** chain/mem/sw (later phases),
`--walk-mode=per-assembly`, and `-L`.

**Code:**
- `pairing.c/h`: pure helpers.
- `search.c`:
  - `refmap_pair_batch` (step 1);
  - `refmap_pair_event` and `write_refmap` (step 2);
  - the pair-aligned reader in `worker_pipeline`.

## Validation (phase 1)

- **Unit:** `test/test_pairing.c` covers mate names, the insert-size estimator (hand-computed),
  every resolver class and tie-break, the near-projection founder filter, and every intersection
  fallback.
- **Integration:** `test/run_integration.sh` section 11 on the `docs/examples/pe/` fixture
  (`test/make_pe_fixture.py` documents its design). It checks:
  - the class of every hand-built pair;
  - P1's Mo17-only event set;
  - P5's rescue past a TE copy;
  - PS4G totals, per-read lines and npy↔PS4G agreement;
  - interleaved vs R1/R2 identity, pair-safe `-K` batches and `--pe-isize`;
  - all input-error paths;
  - a valgrind pass.
- **Real data, single-end identity:** `grits_workdir/scratch/paired_end/se_identity.sh` compares
  stdout, PS4G, npy and sidecars on 100k real reads, for refmap (production placed-fix flags),
  chain, and chain `--support=best --dup-count`, against the parent binaries.
- **Maize wgsim pairs:** `grits_workdir/scratch/paired_end/pe_sim_eval.sh`. It simulates 500k pairs
  each from Oh43 (in-panel) and Ia453 (held-out), runs refmap single-end and `--paired` on the
  same reads, and scores mate placement (`placement_truth.py`) and events
  (`pe_event_eval.py`). Results are below.

## Phase-1 results (2026-10-14, binary `1238c8b`)

**Single-end identity:** with `--paired` off, all 18 outputs are IDENTICAL to the parent binaries.
That is stdout, PS4G, npy and 3 sidecars for each of refmap, chain, and chain
`--support=best --dup-count`, on 100k real reads.

**Maize wgsim pairs:** 500k pairs per sample, 150 bp, inserts 500±50, `-e 0.001`. Production
refmap flags (k50 lift, `--lift-chain-nb=3 --lift-chain-interp`). SE is the same reads
single-end. Placement percentages are per mate row against the truth gVCF map; event
percentages are from `--ps4g-per-read`.

| | Oh43 SE | Oh43 PE | Oh43 PE+trust | Ia453 SE | Ia453 PE | Ia453 PE+trust |
|---|---|---|---|---|---|---|
| EXACT mates, wrong chr | 1.33 | 1.04 | 1.04 | 1.79 | 1.50 | 1.50 |
| PLACED collinear, <2 kb | 76.36 | 76.83 | 77.62 | 73.95 | 74.43 | 75.46 |
| PLACED collinear, wrong chr | 3.03 | 2.50 | 2.23 | 4.48 | 3.88 | 3.45 |
| PLACED non-collinear, wrong chr | 6.97 | 5.38 | 7.02 | 10.41 | 8.79 | 10.52 |
| events per fragment | 1.79 | 1.04 | 1.02 | 1.78 | 1.07 | 1.03 |
| wrong events (≥20 kb or wrong chr), count | 27,267 | 18,597 | 16,368 | 31,391 | 22,790 | 19,625 |
| fragments with any wrong event | 7.17 | 6.28 | 5.62 | 9.26 | 8.34 | 7.37 |
| mean founder-set size | 9.89 | 8.13 | 7.96 | 9.18 | 7.56 | 7.32 |
| true founder in the set | 98.85 | 97.38 | 97.46 | — | — | — |

Pair classes, Oh43 default: 317k CONC, 17k RESC, 98k SOLO (mostly MULTI + one placed mate),
17k DISC, 50k NONE. The estimated insert is mean 495 / sd 52.

**Reading the table:**
- Pairing cuts wrong events by about 32% and makes founder sets about 18% smaller.
- The intersection costs about 1.5 pp of true-founder recall on the in-panel line. In those
  cases one mate's sampled set is missing the true founder.
- 20–100 kb PLACED misplacement does not move. Those reads have no lift projection near
  their mate; the whole projection set is shifted.
- `--pe-trust-exact` is a trade-off. It helps collinear fragments but moves non-collinear
  ones onto EXACT mates at repeat copies.
- **Not yet measured:** the effect on GRITS model accuracy. The model was trained on
  single-end event statistics (see Next phases).

## Next phases

- **Phase 2, chain:**
  - Link SMEM chains across mates: a mate's chain that is colinear with the other's within the
    insert counts as one fragment.
  - One event per pair, with the founder set = intersection.
  - `--dup-count` must dedup on the pair; the `;x=N` suffix would sit on the fragment name.
- **Phase 3, mem/sw:** bwa-style pair scoring and mate rescue (SW of the unaligned mate inside
  the insert window).
- **Downstream:**
  - The GRITS simulator and training corpus are single-end. Paired PS4G events have different
    statistics (fewer events, tighter founder sets), so a model trained on single-end features
    needs a recalibration or retrain check before paired data is used for real.
  - Retune the router `disj` whenever placement changes (memory `placed_fix_chain_projection`).
