#!/usr/bin/env python3
"""Score refmap placements against source coordinates.

Joins refmap output to the ground-truth table (query name encodes
LINE|chr|pos0|len) and classifies each query by whether it anchored back to a
B73 coordinate on the SAME chromosome within a tolerance (default +-5 Mb).

A query's reference position is taken as the breakpoint coordinate(s) cL/cR; a
placement counts as correct when the reference chromosome matches the source
chromosome and min(|cL-pos|,|cR-pos|) <= tol. NAM founders are largely colinear
with B73, so the founder source coordinate is a proxy for the expected B73
coordinate; the tolerance absorbs indel drift.

refmap columns: qname qlen status nCar carriers refName strand cL cR span ins

Usage:
    analyze.py --refmap refmap.out --truth truth.tsv [--tol 5000000]
"""
import argparse, sys
from collections import defaultdict


def chr_of(refname):
    # "B73_chr5" -> "chr5"; "." -> None
    if not refname or refname == ".":
        return None
    return refname.split("_", 1)[1] if "_" in refname else refname


def pct(n, d):
    return f"{100.0*n/d:6.2f}%" if d else "   n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--tol", type=float, default=5e6)
    a = ap.parse_args()
    tol = a.tol

    truth = {}
    with open(a.truth) as f:
        next(f)
        for ln in f:
            q, line, c, pos0, L = ln.rstrip("\n").split("\t")
            truth[q] = (line, c, int(pos0))

    # per-line and per-status tallies
    lines = ["B73", "B97", "Ki3", "CML247"]
    # bucket counts keyed by (line)
    cnt = defaultdict(lambda: defaultdict(int))   # cnt[line][bucket]
    status_mix = defaultdict(lambda: defaultdict(int))  # status_mix[line][status]
    deltas = defaultdict(list)  # deltas[line] -> list of |delta| for same-chr placements
    ncar_sum = defaultdict(int)
    seen = set()
    n_rows = 0

    with open(a.refmap) as f:
        for ln in f:
            if not ln.strip() or ln[0] == "#":
                continue
            F = ln.rstrip("\n").split("\t")
            if len(F) < 11:
                continue
            q, qlen, status, nCar, carriers, refName, strand, cL, cR, span, ins = F[:11]
            if q not in truth:
                continue
            n_rows += 1
            seen.add(q)
            line, src_c, src_pos = truth[q]
            status_mix[line][status] += 1
            try:
                ncar_sum[line] += int(nCar)
            except ValueError:
                pass

            rc = chr_of(refName)
            if status == "UNPLACED" or rc is None or cL == ".":
                cnt[line]["UNPLACED"] += 1
                continue
            if rc != src_c:
                cnt[line]["wrong_chr"] += 1
                continue
            # same chromosome: compute distance from source to nearest breakpoint edge
            edges = []
            for v in (cL, cR):
                if v != ".":
                    edges.append(abs(int(v) - src_pos))
            if not edges:
                cnt[line]["UNPLACED"] += 1
                continue
            d = min(edges)
            deltas[line].append(d)
            if d <= tol:
                cnt[line][f"correct_{status}"] += 1
                cnt[line]["correct"] += 1
            else:
                cnt[line]["right_chr_far"] += 1

    missing = set(truth) - seen
    if missing:
        sys.stderr.write(f"WARNING: {len(missing)} truth queries had no refmap row\n")

    def line_total(line):
        return sum(status_mix[line].values())

    # ---- report ----
    out = []
    out.append(f"# Ref-sensitivity: refmap placement accuracy (tol = +-{tol/1e6:g} Mb)")
    out.append("")
    out.append(f"Total refmap rows joined to truth: {n_rows}")
    if missing:
        out.append(f"Truth queries with no refmap row: {len(missing)}")
    out.append("")
    out.append("## Headline: fraction anchored to correct B73 chr within tolerance")
    out.append("")
    out.append(f"{'line':8} {'n':>7} {'correct':>9} {'%correct':>9} "
               f"{'wrong_chr':>10} {'far(>tol)':>10} {'unplaced':>9}")
    grand = defaultdict(int)
    for line in lines:
        n = line_total(line)
        if n == 0:
            continue
        c = cnt[line]
        out.append(f"{line:8} {n:7d} {c['correct']:9d} {pct(c['correct'],n):>9} "
                   f"{c['wrong_chr']:10d} {c['right_chr_far']:10d} {c['UNPLACED']:9d}")
        for k, v in c.items():
            grand[k] += v
        grand["n"] += n
    out.append(f"{'ALL':8} {grand['n']:7d} {grand['correct']:9d} {pct(grand['correct'],grand['n']):>9} "
               f"{grand['wrong_chr']:10d} {grand['right_chr_far']:10d} {grand['UNPLACED']:9d}")
    out.append("")

    # carriers vs control split
    carriers_lines = [l for l in lines if l != "B73"]
    car_n = sum(line_total(l) for l in carriers_lines)
    car_correct = sum(cnt[l]["correct"] for l in carriers_lines)
    out.append("## Control vs carriers")
    out.append(f"  B73 control : {cnt['B73']['correct']}/{line_total('B73')} "
               f"correct ({pct(cnt['B73']['correct'], line_total('B73'))})")
    out.append(f"  carriers    : {car_correct}/{car_n} correct ({pct(car_correct, car_n)})")
    out.append("")

    out.append("## Status mix per line")
    statuses = ["EXACT", "PLACED", "ONE_SIDE", "UNPLACED"]
    out.append(f"{'line':8} " + " ".join(f"{s:>9}" for s in statuses) + f" {'avgCarr':>8}")
    for line in lines:
        n = line_total(line)
        if n == 0:
            continue
        row = " ".join(f"{status_mix[line].get(s,0):9d}" for s in statuses)
        avgc = ncar_sum[line]/n if n else 0
        out.append(f"{line:8} {row} {avgc:8.2f}")
    out.append("")

    out.append("## Correct-placement breakdown by status")
    out.append(f"{'line':8} {'EXACT_ok':>10} {'PLACED_ok':>10} {'ONE_SIDE_ok':>12}")
    for line in lines:
        if line_total(line) == 0:
            continue
        c = cnt[line]
        out.append(f"{line:8} {c.get('correct_EXACT',0):10d} {c.get('correct_PLACED',0):10d} "
                   f"{c.get('correct_ONE_SIDE',0):12d}")
    out.append("")

    out.append("## Distance distribution |refpos - srcpos| for SAME-chr placements (bp)")
    out.append(f"{'line':8} {'n':>8} {'median':>12} {'p90':>12} {'p99':>12} {'max':>14}")
    def q(sl, frac):
        if not sl:
            return 0
        i = min(len(sl)-1, int(frac*len(sl)))
        return sl[i]
    all_d = []
    for line in lines:
        ds = sorted(deltas[line])
        all_d += ds
        if ds:
            out.append(f"{line:8} {len(ds):8d} {q(ds,0.5):12d} {q(ds,0.9):12d} "
                       f"{q(ds,0.99):12d} {ds[-1]:14d}")
    all_d.sort()
    if all_d:
        out.append(f"{'ALL':8} {len(all_d):8d} {q(all_d,0.5):12d} {q(all_d,0.9):12d} "
                   f"{q(all_d,0.99):12d} {all_d[-1]:14d}")
    out.append("")
    # ECDF at a few thresholds
    out.append("## Same-chr placement within distance thresholds (cumulative, ALL lines)")
    for thr, lbl in [(1e4,"10kb"),(1e5,"100kb"),(1e6,"1Mb"),(5e6,"5Mb"),(1e7,"10Mb"),(5e7,"50Mb")]:
        k = sum(1 for d in all_d if d <= thr)
        out.append(f"  <= {lbl:6}: {k:8d} / {len(all_d)}  ({pct(k,len(all_d))})")

    text = "\n".join(out)
    print(text)


if __name__ == "__main__":
    main()
