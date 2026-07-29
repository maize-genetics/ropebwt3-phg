#!/usr/bin/env python3
"""Calibrate k-mer-agreement confidence signals against measured precision.

For each read: gather its k-mer votes (via the liftover), cluster them, and record
  n_vote   = distinct k-mers that cast any vote (informative tiles)
  best     = distinct k-mers in the winning cluster (agree on the placement)
  second   = distinct k-mers in the runner-up cluster (competing locus)
  margin   = best - second
  frac     = best / n_vote
then join truth and measure precision within bins of each signal. If a signal is
well calibrated, precision rises monotonically with it -> it can be reported as a
per-read confidence / MAPQ.
"""
import argparse, sys, statistics
from collections import defaultdict
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from liftover import parse_pos, load_anchors, make_index, project_robust


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors-mem", required=True)
    ap.add_argument("--kmers-mem", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--k", type=int, default=75)
    ap.add_argument("--max-occ", type=int, default=4)
    ap.add_argument("--cluster-tol", type=float, default=2000)
    ap.add_argument("--tol", type=float, default=5e6)
    ap.add_argument("--ref", default="B73")
    a = ap.parse_args()

    truth = {}
    with open(a.truth) as f:
        next(f)
        for ln in f:
            q, line, c, p, L = ln.rstrip().split("\t")
            truth[q] = (line, c, int(p))
    lift = load_anchors(a.anchors_mem, a.max_occ, a.ref)
    idx = make_index(lift)

    votes = defaultdict(list)  # read -> [(chr,pos,kid)]
    with open(a.kmers_mem) as f:
        for ln in f:
            F = ln.rstrip("\n").split("\t")
            if len(F) < 5:
                continue
            read, off = F[0].rsplit("#", 1)
            if read not in truth or int(F[2]) - int(F[1]) < a.k or int(F[3]) > a.max_occ:
                continue
            got = set()
            for t in F[5:]:
                tx, c, st, p = parse_pos(t)
                if tx == a.ref:
                    got.add((c, p))
                else:
                    r = project_robust(idx, tx, c, p, 500000, max_mad=200000)
                    if r is not None:
                        got.add(r)
            for (c, p) in got:
                votes[read].append((c, p, off))

    # per-read stats
    rows = []  # (n_vote, best, second, correct_or_None)
    for read, (line, sc, sp) in truth.items():
        vs = votes.get(read, [])
        n_vote = len(set(k for _, _, k in vs))
        # cluster per chr
        clusters = []  # (support, chr, medpos)
        byc = defaultdict(list)
        for (c, p, k) in vs:
            byc[c].append((p, k))
        for c, lst in byc.items():
            lst.sort()
            i = 0
            while i < len(lst):
                j = i; kids = set(); poss = []
                while j < len(lst) and lst[j][0] - lst[i][0] <= a.cluster_tol:
                    kids.add(lst[j][1]); poss.append(lst[j][0]); j += 1
                clusters.append((len(kids), c, int(statistics.median(poss))))
                i = j
        clusters.sort(reverse=True)
        best = clusters[0][0] if clusters else 0
        second = clusters[1][0] if len(clusters) > 1 else 0
        if best == 0:
            rows.append((n_vote, best, second, None)); continue
        pc, pchr, ppos = clusters[0]
        ok = (pchr == sc and abs(ppos - sp) <= a.tol)
        rows.append((n_vote, best, second, ok))

    def report(name, keyfn):
        print(f"\n## precision by {name}")
        b = defaultdict(lambda: [0, 0])  # key -> [placed, correct]
        for n_vote, best, second, ok in rows:
            if ok is None:
                continue
            k = keyfn(n_vote, best, second)
            b[k][0] += 1; b[k][1] += int(ok)
        print(f"  {name:22} {'placed':>8} {'precision':>10}")
        for k in sorted(b):
            pl, co = b[k]
            print(f"  {str(k):22} {pl:8d} {100.0*co/pl:9.1f}%")

    placed = sum(1 for r in rows if r[3] is not None)
    correct = sum(1 for r in rows if r[3])
    print(f"reads={len(rows)} placed(best>=1)={placed} overall_prec={100.0*correct/placed:.1f}%")
    report("best support (agree)", lambda nv, b, s: b)
    report("margin best-second", lambda nv, b, s: b - s)
    report("agreement frac best/n_vote", lambda nv, b, s: f"{b}/{nv}" if nv <= 6 else ">6")


if __name__ == "__main__":
    main()
