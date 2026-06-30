#!/usr/bin/env python3
"""E4 prototype, step 2: build a carrier->B73 liftover from located anchors and
self-validate it.

Input: `mem -p` output on the anchors (name = CARRIER|chr|pos|len). For every
anchor that is single-copy-per-taxon (count <= max_occ) and has exactly one B73
hit, we get a pair (carrier_chr, carrier_pos) <-> (B73_chr, B73_pos). These are
sorted per (carrier, carrier_chr) into a piecewise-linear liftover.

project(carrier, chr, pos): find the flanking anchors; if both map to the SAME
B73 chr and are locally collinear (monotone, slope within tol of 1), linearly
interpolate -> approximate B73 pos. Otherwise return None (NULL slot) -- we would
rather not place than place wrong.

Self-test: leave-one-out over the anchors themselves -- project each anchor from
its neighbours and compare to its own known B73 position.
"""
import argparse, sys, bisect
from collections import defaultdict


def parse_pos(tok):
    # "B73_chr2:+:40493267" -> (taxon, chr, strand, pos)
    name, strand, pos = tok.rsplit(":", 2)
    taxon = name.split("_", 1)[0]
    chrom = name.split("_", 1)[1] if "_" in name else name
    return taxon, chrom, strand, int(pos)


def load_anchors(path, max_occ, ref="B73"):
    # per (carrier, carrier_chr): list of (carrier_pos, b73_chr, b73_pos)
    lift = defaultdict(list)
    n_total = n_used = 0
    with open(path) as f:
        for ln in f:
            F = ln.rstrip("\n").split("\t")
            if len(F) < 5:
                continue
            n_total += 1
            qname = F[0]
            count = int(F[3])
            if count > max_occ:
                continue
            car, cchr, cpos, _ = qname.split("|")
            cpos = int(cpos)
            toks = F[5:]
            b73 = [parse_pos(t) for t in toks if t.startswith(ref + "_")]
            if len(b73) != 1:
                continue
            _, bchr, bstr, bpos = b73[0]
            lift[(car, cchr)].append((cpos, bchr, bpos))
            n_used += 1
    for k in lift:
        lift[k].sort()
    sys.stderr.write(f"anchors: {n_total} total, {n_used} usable liftover points "
                     f"({100.0*n_used/n_total:.1f}%)\n")
    return lift


def make_index(lift):
    # for bisect: per key, parallel arrays
    idx = {}
    for k, arr in lift.items():
        cpos = [a[0] for a in arr]
        idx[k] = (cpos, arr)
    return idx


def project(idx, car, cchr, pos, slope_tol=0.5, max_gap=2_000_000):
    key = (car, cchr)
    if key not in idx:
        return None
    cpos, arr = idx[key]
    j = bisect.bisect_left(cpos, pos)
    # flanking anchors lo (<=pos) and hi (>pos)
    lo = j - 1
    hi = j if j < len(arr) and arr[j][0] != pos else j
    if j < len(arr) and arr[j][0] == pos:
        # exact anchor hit
        return arr[j][1], arr[j][2]
    if lo < 0 or hi >= len(arr):
        return None
    p0, c0, b0 = arr[lo]
    p1, c1, b1 = arr[hi]
    if c0 != c1:
        return None                      # flanks disagree on B73 chr -> NULL
    if p1 - p0 == 0 or p1 - p0 > max_gap:
        return None                      # too sparse to trust
    slope = (b1 - b0) / (p1 - p0)
    if abs(abs(slope) - 1.0) > slope_tol:
        return None                      # not locally collinear -> NULL
    bpos = b0 + slope * (pos - p0)
    return c1, int(bpos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors-mem", required=True)
    ap.add_argument("--max-occ", type=int, default=4)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    lift = load_anchors(a.anchors_mem, a.max_occ)
    idx = make_index(lift)
    # density report
    keys = sorted(lift.keys())
    npts = sum(len(v) for v in lift.values())
    sys.stderr.write(f"liftover keys (carrier,chr): {len(keys)}; points: {npts}\n")
    if a.selftest:
        # leave-one-out: project each anchor from its neighbours
        import statistics
        ok = nulls = far = n = 0
        deltas = []
        for (car, cchr), arr in lift.items():
            cpos = [x[0] for x in arr]
            for i in range(1, len(arr) - 1):
                pos, bchr, bpos = arr[i]
                # temporarily project from neighbours i-1, i+1
                p0, c0, b0 = arr[i - 1]
                p1, c1, b1 = arr[i + 1]
                n += 1
                if c0 != c1 or c0 != bchr:
                    nulls += 1
                    continue
                if p1 == p0:
                    nulls += 1; continue
                slope = (b1 - b0) / (p1 - p0)
                if abs(abs(slope) - 1.0) > 0.5:
                    nulls += 1; continue
                est = b0 + slope * (pos - p0)
                d = abs(est - bpos)
                deltas.append(d)
                if d <= 5e6:
                    ok += 1
                else:
                    far += 1
        deltas.sort()
        def q(f): return int(deltas[min(len(deltas)-1, int(f*len(deltas)))]) if deltas else 0
        sys.stderr.write(
            f"SELFTEST leave-one-out: n={n} projected={len(deltas)} null={nulls} "
            f"({100.0*nulls/n:.1f}%) within5Mb={ok} far={far}\n")
        if deltas:
            sys.stderr.write(
                f"  |delta B73| median={q(0.5)} p90={q(0.9)} p99={q(0.99)} max={deltas[-1]}\n")


if __name__ == "__main__":
    main()
