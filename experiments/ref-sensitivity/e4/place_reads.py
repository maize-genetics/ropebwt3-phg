#!/usr/bin/env python3
"""E4 prototype, step 3: place reads via the carrier->B73 liftover and score.

For each read (from `mem -p` output):
  occ > max_occ                         -> MULTI    (repeat/retro, not placed)
  has a B73 hit (occ <= max_occ)        -> EXACT    (place at the B73 coordinate)
  carrier-only                          -> project each carrier hit via the
                                           liftover; if >=1 projects and they
                                           agree on chr (spread <= agree_tol)
                                           -> PLACED at the median; else NULL/UNPLACED

Scored against truth.tsv (source coord as B73 proxy, tol +-5 Mb), reporting the
same precision/recall/wrong-chr as the refmap E0-E3 runs for direct comparison.
"""
import argparse, sys, bisect, statistics
from collections import defaultdict
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from liftover import parse_pos, load_anchors, make_index, project, project_robust


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors-mem", required=True)
    ap.add_argument("--reads-mem", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--max-occ", type=int, default=4)
    ap.add_argument("--agree-tol", type=float, default=5e6)
    ap.add_argument("--tol", type=float, default=5e6)
    ap.add_argument("--ref", default="B73")
    ap.add_argument("--robust", action="store_true", help="windowed slope=+/-1 robust projection")
    ap.add_argument("--win", type=float, default=500000)
    ap.add_argument("--max-mad", type=float, default=200000)
    a = ap.parse_args()
    proj_fn = (lambda car, c, p: project_robust(idx, car, c, p, int(a.win), max_mad=a.max_mad)) if a.robust \
              else (lambda car, c, p: project(idx, car, c, p))

    truth = {}
    with open(a.truth) as f:
        next(f)
        for ln in f:
            q, line, c, p, L = ln.rstrip().split("\t")
            truth[q] = (line, c, int(p))

    lift = load_anchors(a.anchors_mem, a.max_occ, a.ref)
    idx = make_index(lift)

    smix = defaultdict(int)
    n = correct = wrong = placed = 0
    b73n = b73ok = 0
    cls = defaultdict(lambda: [0, 0, 0])  # status -> [n, correct, wrong]

    with open(a.reads_mem) as f:
        for ln in f:
            F = ln.rstrip("\n").split("\t")
            if len(F) < 5:
                continue
            q = F[0]
            if q not in truth:
                continue
            line, sc, sp = truth[q]
            n += 1
            isb = (line == a.ref)
            if isb:
                b73n += 1
            count = int(F[3])
            toks = F[5:]
            poss = [parse_pos(t) for t in toks]
            status = pchr = None
            ppos = -1
            if count > a.max_occ:
                status = "MULTI"
            else:
                b73 = [(c, p) for (tx, c, st, p) in poss if tx == a.ref]
                if b73:
                    chrs = set(c for c, _ in b73)
                    if len(chrs) == 1:
                        status, pchr, ppos = "EXACT", b73[0][0], b73[0][1]
                    else:
                        status = "AMBIG"     # multiple B73 chrs: untrusted
                else:
                    # carrier-only: project each carrier hit
                    proj = []
                    for (tx, c, st, p) in poss:
                        if tx == a.ref:
                            continue
                        r = proj_fn(tx, c, p)
                        if r is not None:
                            proj.append(r)
                    if proj:
                        pc = defaultdict(list)
                        for c, p in proj:
                            pc[c].append(p)
                        bestc = max(pc, key=lambda c: len(pc[c]))
                        ps = sorted(pc[bestc])
                        if ps[-1] - ps[0] <= a.agree_tol:
                            status, pchr, ppos = "PLACED", bestc, int(statistics.median(ps))
                        else:
                            status = "UNPLACED"   # carrier projections disagree
                    else:
                        status = "UNPLACED"       # all NULL
            smix[status] += 1
            if status in ("EXACT", "PLACED"):
                placed += 1
                ok = (pchr == sc and abs(ppos - sp) <= a.tol)
                if pchr != sc:
                    wrong += 1
                if ok:
                    correct += 1
                    if isb:
                        b73ok += 1
                t = cls[status]; t[0] += 1; t[1] += int(ok); t[2] += int(pchr != sc)

    prec = 100.0 * correct / placed if placed else 0
    rec = 100.0 * correct / n if n else 0
    print(f"reads={n}  placed={placed}  correct={correct}  wrong_chr={wrong}")
    print(f"PRECISION={prec:.1f}%  RECALL={rec:.1f}%  B73ctrl={100.0*b73ok/b73n:.1f}%")
    print("status mix: " + "  ".join(f"{k}={v}" for k, v in sorted(smix.items())))
    print("by-status (n / correct / wrong_chr):")
    for k in ("EXACT", "PLACED"):
        t = cls[k]
        if t[0]:
            print(f"  {k:8} n={t[0]:6d} correct={t[1]:6d} ({100.0*t[1]/t[0]:.1f}%) wrong_chr={t[2]}")


if __name__ == "__main__":
    main()
