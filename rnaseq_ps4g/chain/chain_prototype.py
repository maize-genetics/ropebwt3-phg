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
        gametes, refs = set(), []
        for t in toks:
            name, strand, pos = t.rsplit(":", 2)
            gametes.add(name2idx[simlib.gamete_of(name)])
            if name.startswith(ref_prefix):
                refs.append((int(pos), strand))   # all reference occurrences (tandem dup -> >1)
        smem = {"qs": qs, "qe": qe, "count": count, "gametes": gametes,
                "refs": refs, "ref": refs[0] if refs else None,
                "informative": count <= max_occ}
        if rid != cur:
            if cur is not None:
                yield cur, rows
            cur, rows = rid, []
        rows.append(smem)
    if cur is not None:
        yield cur, rows


def canonical_splice(refseq, strand, lo, hi, slack=6):
    """Does the reference intron near [lo,hi) have canonical splice motifs (GT..AG
    in transcription orientation; '-' strand checked on the reverse complement)?

    The exon/intron boundary implied by a SMEM's end drifts a few bp from the true
    splice site (microhomology: an intron base coincidentally matches the adjacent
    exon, so the SMEM over/under-runs). Like real spliced aligners we therefore
    search a small ±`slack` window on each boundary for the canonical motif. A
    non-canonical intron-sized gap is a chimera/artifact, not a real splice."""
    if hi - lo < 4:
        return True                       # too small to be a real intron -> allow
    n = len(refseq)
    # transcription-orientation GT..AG on the forward strand reads GT..AG for a '+'
    # intron and CT..AC (revcomp) for a '-' intron.
    low_motif, high_motif = ("GT", "AG") if strand == "+" else ("CT", "AC")
    for d in range(max(0, lo - slack), min(n - 1, lo + slack) + 1):
        if refseq[d:d + 2] != low_motif:
            continue
        for a in range(max(d + 4, hi - slack), min(n, hi + slack) + 1):
            if refseq[a - 2:a] == high_motif:
                return True
    return False


def chain_and_intersect(smems, gap_intron, max_intron=500, gap_coef=0.02,
                        refseq=None, splice_min=10):
    """Return (segments, gset, flag) for one read. flag is 'ok', 'unplaced' (no
    chain) or 'ambiguous' (maps to >1 reference locus -- not confidently placeable).

    - Keep informative SMEMs with a reference hit; pick the strand and the maximal
      IN-ORDER colinear chain by a `(#anchors, bases - gap_penalty)` DP.
    - The gap penalty is on the *unexplained* reference jump `Δref - Δquery` (0 for
      a contiguous within-exon link; = the intron length for a splice). It grows
      with the jump, so among equal-anchor chains a compact one beats a distant
      (chimeric/paralogous) link, while a real large intron still wins when it is
      the only multi-anchor option (anchors dominate the score). Links whose
      unexplained jump exceeds `max_intron` are rejected outright.
    - Intersect the gamete sets over the whole chain (junction linkage), and split
      into exon segments at reference gaps > gap_intron.
    - Ambiguity: if a chained SMEM has >1 reference occurrence (e.g. a tandem
      duplication), the read maps to multiple loci -> flag 'ambiguous', don't emit
      a confident position (suppressing false evidence is safer than guessing).
    """
    cand = [s for s in smems if s["informative"] and s["ref"] is not None]
    if not cand:
        return None, None, "unplaced"
    # strand = the one carrying the most query bases (ties -> '+')
    span = {"+": 0, "-": 0}
    for s in cand:
        span[s["ref"][1]] += s["qe"] - s["qs"]
    strand = "+" if span["+"] >= span["-"] else "-"
    cand = [s for s in cand if s["ref"][1] == strand]
    cand.sort(key=lambda s: (s["qs"], s["ref"][0]))
    n = len(cand)
    w = [s["qe"] - s["qs"] for s in cand]
    # best[i] = (n_anchors, bases - gap_penalty) of the best chain ending at i
    best = [(1, float(w[i])) for i in range(n)]
    prev = [-1] * n
    for i in range(n):
        ri = cand[i]["ref"][0]
        for j in range(i):
            rj = cand[j]["ref"][0]
            dr = (ri - rj) if strand == "+" else (rj - ri)   # forward ref advance
            if dr < 0:                                        # out of order -> not colinear
                continue
            dq = max(0, cand[i]["qs"] - cand[j]["qe"])        # query advance (clamp overlap)
            unexplained = max(0, dr - dq)                     # intron length / paralog jump
            if unexplained > max_intron:                      # implausible -> reject the link
                continue
            # GT-AG splice-site check: an intron-sized gap must have canonical
            # boundaries, else it is a non-canonical (chimeric) junction -> reject.
            if refseq is not None:
                wj, wi = cand[j]["qe"] - cand[j]["qs"], cand[i]["qe"] - cand[i]["qs"]
                lo, hi = (rj + wj, ri) if strand == "+" else (ri + wi, rj)
                if hi - lo >= splice_min and not canonical_splice(refseq, strand, lo, hi):
                    continue
            score = (best[j][0] + 1, best[j][1] + w[i] - gap_coef * unexplained)
            if cand[j]["qs"] <= cand[i]["qs"] and score > best[i]:
                best[i] = score
                prev[i] = j
    end = max(range(n), key=lambda i: best[i])
    idx, chain = end, []
    while idx != -1:
        chain.append(cand[idx])
        idx = prev[idx]
    chain.reverse()
    if not chain:
        return None, None, "unplaced"
    if any(len(s["refs"]) > 1 for s in chain):                # >1 reference locus
        return None, None, "ambiguous"
    gset = set.intersection(*[s["gametes"] for s in chain])
    ivs = sorted((s["ref"][0], s["ref"][0] + (s["qe"] - s["qs"])) for s in chain)
    segs, lo, hi = [], ivs[0][0], ivs[0][1]
    for a, b in ivs[1:]:
        if a - hi > gap_intron:               # intron-sized gap -> new exon segment
            segs.append((lo, hi))
            lo, hi = a, b
        else:
            hi = max(hi, b)
    segs.append((lo, hi))
    return segs, gset, "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gametes", required=True)
    ap.add_argument("--ref-prefix", default="B73")
    ap.add_argument("--contig", default="chr1")
    ap.add_argument("--max-occ", type=int, default=5, help="SMEMs with count > this are repeats (skipped)")
    ap.add_argument("--gap-intron", type=int, default=30, help="reference gap starting a new exon segment")
    ap.add_argument("--max-intron", type=int, default=500,
                    help="reject a chain link whose unexplained reference jump "
                         "(Δref - Δquery) exceeds this (plausible intron ceiling). "
                         "Conservative default trades junction power for specificity; "
                         "raise for large-intron organisms")
    ap.add_argument("--gap-coef", type=float, default=0.02,
                    help="penalty per bp of unexplained reference jump; disfavors "
                         "distant/chimeric links vs a compact chain of equal anchor count")
    ap.add_argument("--ref-fasta",
                    help="reference (or pangenome) FASTA; enables the GT-AG splice-"
                         "site check so an intron-gap link must have canonical motifs")
    ap.add_argument("--splice-min", type=int, default=10,
                    help="reference gap size above which the GT-AG check applies")
    a = ap.parse_args()

    name2idx = parse_gametes(a.gametes)
    refseq = None
    if a.ref_fasta:
        fa = simlib.read_fasta(a.ref_fasta)
        refseq = next((v for k, v in fa.items() if k.startswith(a.ref_prefix)), None)
    sys.stdout.write("readName\trefContig\trefPos\tgameteSet\n")
    n_reads = n_emit = n_unplaced = n_ambig = 0
    for rid, smems in read_smems(sys.stdin, a.ref_prefix, name2idx, a.max_occ):
        n_reads += 1
        segs, gset, flag = chain_and_intersect(smems, a.gap_intron, a.max_intron,
                                               a.gap_coef, refseq, a.splice_min)
        if flag == "ambiguous":
            n_ambig += 1                      # maps to >1 locus -> suppress (no false evidence)
            continue
        if segs is None or not gset:
            n_unplaced += 1
            continue
        gs = ",".join(map(str, sorted(gset)))
        for lo, hi in segs:
            sys.stdout.write("%s\t%s\t%d\t%s\n" % (rid, a.contig, lo, gs))
            n_emit += 1
    sys.stderr.write("reads=%d emitted_rows=%d unplaced=%d ambiguous=%d\n"
                     % (n_reads, n_emit, n_unplaced, n_ambig))


if __name__ == "__main__":
    main()
