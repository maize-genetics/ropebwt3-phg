#!/usr/bin/env python3
"""Prototype: unite a read's SMEMs by colinear chaining + strict set intersection.

Reads `ropebwt3 mem -p N` output and emits a per-read PS4G file
(readName, refContig, refPos, gameteSet), one row per exon SEGMENT, carrying the
whole-read intersection of the informative colinear SMEM gamete sets.

Why intersection is safe (see the plan): each SMEM is an error-free substring of
the read, which came from the source, so the source is in every SMEM's set and
always survives the intersection. Off-target/chimeric SMEMs are not colinear with
the reference chain and are dropped BEFORE intersecting; repeat SMEMs (count >
max_occ) carry no trustworthy set and are skipped. So intersection only removes
spurious founders (consistent with some SMEMs but not all) -> founder-dropout
stays 0.

Stage-0 scope: reads come from the reference (B73), so every true exon SMEM has a
reference hit; we chain on those. Carrier-only SMEMs (insertions/PAV, or chimeric
junction-spanning matches like the r000001 Mo17 hit) have no reference hit and are
excluded here -- lift projection for genuine carrier-only SMEMs is a later step.

Usage:
  ropebwt3 mem -l 19 -p 16 idx.fmd reads.fq | \
    chain_prototype.py --gametes gametes.tsv --ref-prefix B73 > reads.ps4g
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sim"))
import simlib


def parse_gametes(path):
    m = {}
    for line in open(path):
        idx, name = line.rstrip("\n").split("\t")
        m[name] = int(idx)
    return m


def read_smems(fp, ref_prefix, name2idx, max_occ):
    """Yield (read_name, [smem,...]). Each smem is a dict with the query span, the
    reference (B73) forward-strand position + strand if present, its gamete set,
    and whether it is informative (count <= max_occ)."""
    cur, rows = None, []
    for line in fp:
        F = line.rstrip("\n").split("\t")
        if len(F) < 4:
            continue
        rid, qs, qe, count = F[0], int(F[1]), int(F[2]), int(F[3])
        toks = F[5:] if len(F) > 5 else []
        gametes, ref = set(), None
        for t in toks:
            name, strand, pos = t.rsplit(":", 2)
            gametes.add(name2idx[simlib.gamete_of(name)])
            if name.startswith(ref_prefix):
                ref = (int(pos), strand)          # reference forward-strand coord
        smem = {"qs": qs, "qe": qe, "count": count, "gametes": gametes,
                "ref": ref, "informative": count <= max_occ}
        if rid != cur:
            if cur is not None:
                yield cur, rows
            cur, rows = rid, []
        rows.append(smem)
    if cur is not None:
        yield cur, rows


def chain_and_intersect(smems, gap_intron, max_ref_span=2000):
    """Return (segments, gset) for one read, or (None, None) if unplaceable.
    - Keep informative SMEMs with a reference hit; pick the strand and the maximal
      colinear (monotonic query<->ref) subset.
    - Intersect their gamete sets over the WHOLE read (junction linkage: a founder
      must be exact across every colinear exon of this read).
    - Split the chain into exon segments at reference gaps > gap_intron.
    Each segment is (ref_lo, ref_hi); emit the intersection set at each.
    """
    cand = [s for s in smems if s["informative"] and s["ref"] is not None]
    if not cand:
        return None, None
    # strand = the one carrying the most query bases (ties -> '+')
    span = {"+": 0, "-": 0}
    for s in cand:
        span[s["ref"][1]] += s["qe"] - s["qs"]
    strand = "+" if span["+"] >= span["-"] else "-"
    cand = [s for s in cand if s["ref"][1] == strand]
    # colinear + IN ORDER: sort by query start, then keep the maximal chain whose
    # reference position moves monotonically with query (up for '+', down for '-').
    # DP (longest chain weighted by covered query bases) so a single out-of-order
    # SMEM is dropped rather than allowed to derail the chain.
    cand.sort(key=lambda s: (s["qs"], s["ref"][0]))
    n = len(cand)
    w = [s["qe"] - s["qs"] for s in cand]
    # score a chain by (#anchors, covered bases): a colinear MULTI-anchor chain
    # (agreement) beats a lone long paralog/chimera, matching --kmer's min-agree
    # idea. best[i] = (anchors, bases) of the best chain ending at i.
    best = [(1, w[i]) for i in range(n)]
    prev = [-1] * n
    for i in range(n):
        ri = cand[i]["ref"][0]
        for j in range(i):
            rj = cand[j]["ref"][0]
            ordered = (ri >= rj) if strand == "+" else (ri <= rj)  # ref monotonic with query
            # spatially compact: a colinear read spans ~read-length (+ introns), not
            # a distant locus -> reject a far-away (chimeric) SMEM even if monotonic.
            compact = abs(ri - rj) <= max_ref_span
            cand_score = (best[j][0] + 1, best[j][1] + w[i])
            if cand[j]["qs"] <= cand[i]["qs"] and ordered and compact \
                    and cand_score > best[i]:
                best[i] = cand_score
                prev[i] = j
    end = max(range(n), key=lambda i: best[i])
    idx, chain = end, []
    while idx != -1:
        chain.append(cand[idx])
        idx = prev[idx]
    chain.reverse()
    if not chain:
        return None, None
    gset = set.intersection(*[s["gametes"] for s in chain])
    # reference intervals (forward-strand) of the chained SMEMs, then segment
    ivs = sorted((s["ref"][0], s["ref"][0] + (s["qe"] - s["qs"])) for s in chain)
    segs, lo, hi = [], ivs[0][0], ivs[0][1]
    for a, b in ivs[1:]:
        if a - hi > gap_intron:               # intron-sized gap -> new exon segment
            segs.append((lo, hi))
            lo, hi = a, b
        else:
            hi = max(hi, b)
    segs.append((lo, hi))
    return segs, gset


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gametes", required=True)
    ap.add_argument("--ref-prefix", default="B73")
    ap.add_argument("--contig", default="chr1")
    ap.add_argument("--max-occ", type=int, default=5, help="SMEMs with count > this are repeats (skipped)")
    ap.add_argument("--gap-intron", type=int, default=30, help="reference gap starting a new exon segment")
    ap.add_argument("--max-ref-span", type=int, default=2000,
                    help="max reference distance between chained SMEMs (rejects "
                         "distant/chimeric hits; raise for large introns)")
    a = ap.parse_args()

    name2idx = parse_gametes(a.gametes)
    sys.stdout.write("readName\trefContig\trefPos\tgameteSet\n")
    n_reads = n_emit = n_unplaced = 0
    for rid, smems in read_smems(sys.stdin, a.ref_prefix, name2idx, a.max_occ):
        n_reads += 1
        segs, gset = chain_and_intersect(smems, a.gap_intron, a.max_ref_span)
        if segs is None or not gset:
            n_unplaced += 1
            continue
        gs = ",".join(map(str, sorted(gset)))
        for lo, hi in segs:
            sys.stdout.write("%s\t%s\t%d\t%s\n" % (rid, a.contig, lo, gs))
            n_emit += 1
    sys.stderr.write("reads=%d emitted_rows=%d unplaced=%d\n" % (n_reads, n_emit, n_unplaced))


if __name__ == "__main__":
    main()
