#!/usr/bin/env python3
"""Place a read from its 75 bp k-mers by AGREEMENT.

Each k-mer of a read is located (mem -p) and projected to the reference (EXACT if
it hits the reference, else via the carrier->reference liftover). A read is placed
only where >= min_agree distinct k-mers vote for the same reference locus (cluster
within cluster_tol on one chromosome) -- a sub-read concordance filter. k-mers
that are repeats (occ > max_occ) are ignored; reads with no agreeing cluster are
UNPLACED. This tolerates sequencing error (a k-mer can be error-free when the read
is not) and raises precision (a lone paralogous k-mer is outvoted).

Scored against truth.tsv (read source coord, tol +-5 Mb).
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
    ap.add_argument("--min-agree", type=int, default=2)
    ap.add_argument("--cluster-tol", type=float, default=5000)
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

    # gather per-read votes: read -> list of (chr, pos, kmer_id)
    votes = defaultdict(list)
    reads_seen = set()
    with open(a.kmers_mem) as f:
        for ln in f:
            F = ln.rstrip("\n").split("\t")
            if len(F) < 5:
                continue
            kname = F[0]
            read, off = kname.rsplit("#", 1)
            reads_seen.add(read)
            if read not in truth:
                continue
            if int(F[2]) - int(F[1]) < a.k:      # only full-length exact k-mer hits;
                continue                          # an error-containing k-mer gives partial SMEMs -> skip
            count = int(F[3])
            if count > a.max_occ:
                continue
            kid = off
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
                votes[read].append((c, p, kid))

    n = correct = wrong = placed = b73n = b73ok = unplaced = 0
    for read, (line, sc, sp) in truth.items():
        n += 1
        isb = (line == a.ref)
        if isb:
            b73n += 1
        vs = votes.get(read, [])
        # cluster votes per chromosome
        best_support = 0; best_pos = None; best_chr = None
        byc = defaultdict(list)
        for (c, p, kid) in vs:
            byc[c].append((p, kid))
        for c, lst in byc.items():
            lst.sort()
            i = 0
            while i < len(lst):
                j = i; kids = set(); poss = []
                while j < len(lst) and lst[j][0] - lst[i][0] <= a.cluster_tol:
                    kids.add(lst[j][1]); poss.append(lst[j][0]); j += 1
                if len(kids) > best_support:
                    best_support = len(kids); best_pos = int(statistics.median(poss)); best_chr = c
                i = j
        if best_support >= a.min_agree:
            placed += 1
            if best_chr != sc:
                wrong += 1
            elif abs(best_pos - sp) <= a.tol:
                correct += 1
                if isb:
                    b73ok += 1
        else:
            unplaced += 1

    prec = 100.0 * correct / placed if placed else 0
    rec = 100.0 * correct / n if n else 0
    print(f"reads={n} placed={placed} correct={correct} wrong_chr={wrong} unplaced={unplaced}")
    print(f"PRECISION={prec:.1f}%  RECALL={rec:.1f}%  B73ctrl={100.0*b73ok/b73n if b73n else 0:.1f}%")


if __name__ == "__main__":
    main()
