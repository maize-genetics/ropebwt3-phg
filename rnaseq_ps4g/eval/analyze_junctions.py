#!/usr/bin/env python3
"""Characterize how an aligner handles exon-spanning (junction) reads.

A spliced read covers N exon segments in the reference. Per-segment PS4G emission
(HANDOFF 2.2) wants evidence at *every* segment; stock refmap emits one placement
per read (its longest exact core), so only one segment of a junction read gets
evidence and the rest of that read's exonic support is lost. This script measures
that gap from the per-read PS4G file (refmap --ps4g-per-read).

Usage:
  analyze_junctions.py --truth truth.tsv --ps4g-per-read reads.ps4g --gametes gametes.tsv
"""
import argparse


def parse_gametes(path):
    m = {}
    for line in open(path):
        idx, name = line.rstrip("\n").split("\t")
        m[name] = int(idx)
    return m


def parse_truth(path):
    out = {}
    with open(path) as f:
        f.readline()
        for line in f:
            F = line.rstrip("\n").split("\t")
            rid, src, njunc, segstr = F[0], F[1], int(F[5]), F[4]
            segs = []
            for seg in segstr.split(";"):
                contig, span = seg.rsplit(":", 1)
                lo, hi = span.split("-")
                segs.append((contig, int(lo), int(hi)))
            out[rid] = {"src": src, "njunc": njunc, "segs": segs}
    return out


def parse_side(path):
    """readName -> list of (contig, pos, set). Stock refmap emits one row per read;
    the chaining emitter emits one per exon segment."""
    out = {}
    with open(path) as f:
        f.readline()
        for line in f:
            rid, contig, pos, gs = line.rstrip("\n").split("\t")
            out.setdefault(rid, []).append(
                (contig, int(pos), {int(x) for x in gs.split(",") if gs}))
    return out


def seg_of(pos, contig, segs):
    """Index of the true segment containing the emitted base, else -1."""
    for i, (c, lo, hi) in enumerate(segs):
        if c == contig and lo <= pos < hi:
            return i
    return -1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--ps4g-per-read", dest="ps4g_per_read", required=True)
    ap.add_argument("--gametes", required=True)
    a = ap.parse_args()

    name2idx = parse_gametes(a.gametes)
    truth = parse_truth(a.truth)
    side = parse_side(a.ps4g_per_read)

    # bucket by number of junctions the read spans
    classes = [("within-exon (0 junc)", lambda j: j == 0),
               ("1 junction", lambda j: j == 1),
               ("2+ junctions", lambda j: j >= 2)]
    for label, pred in classes:
        rids = [r for r, t in truth.items() if pred(t["njunc"])]
        if not rids:
            continue
        n = len(rids)
        seg_instances = 0     # total true exon segments across these reads
        near = 0              # read placed near a true exon (some row's [pos,pos+L) overlaps a segment)
        covered = 0           # true exon segments with an emitted base inside them (valid-base coverage)
        recall_iv = 0         # placed near a true exon AND source founder in the emitted set
        for rid in rids:
            t = truth[rid]
            segs = t["segs"]
            seg_instances += len(segs)
            L = sum(hi - lo for _, lo, hi in segs) or 1   # read length in ref bases
            rows = side.get(rid)
            if not rows:
                continue
            allset = set().union(*[r[2] for r in rows])
            # a segment is covered if some emitted base lands inside it
            for (c, lo, hi) in segs:
                if any(rc == c and lo <= rp < hi for (rc, rp, _) in rows):
                    covered += 1
            if any(rc == c and rp < hi and lo < rp + L
                   for (rc, rp, _) in rows for (c, lo, hi) in segs):
                near += 1
                if name2idx.get(t["src"]) in allset:
                    recall_iv += 1
        pct = lambda x, d: 100.0 * x / d if d else 0.0
        print("%-22s n=%-5d segments/read=%.2f  recall(near-exon)=%.1f%%" % (
            label, n, seg_instances / n, pct(recall_iv, n)))
        print("%-22s   segment coverage (valid emitted base) = %.1f%% (%d/%d segments)" % (
              "", pct(covered, seg_instances), covered, seg_instances))
    print("\nSegment coverage = fraction of true exon segments that receive an emitted")
    print("base. Stock refmap emits one position/read (capped at 1/N for a spliced read,")
    print("and its extrapolated base often lands off the exons); the chaining emitter emits")
    print("each colinear SMEM at its own exon -> covers all N exons at valid bases.")


if __name__ == "__main__":
    main()
