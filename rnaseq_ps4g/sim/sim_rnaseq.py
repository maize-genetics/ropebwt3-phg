#!/usr/bin/env python3
"""Module 0, part 2: transcribe/splice the synthetic pangenome into RNAseq reads
with exact ground truth.

Consumes sim_genomes.py output (sim.json + pangenome.fa). For each source
genotype it builds mature mRNAs (exons spliced in transcription order -> real
exon-exon junctions), draws transcripts from a log-normal expression
distribution, fragments them into fixed-length reads, and applies a
substitution error model. Truth is exact by construction because every read
position maps back through the genotype's EditMap to a reference coordinate.

Stage-0 defaults (HANDOFF milestone 0): source = reference only, substitution
error only, single-end, sense strand, within-exon reads (--max-junctions 0).

Outputs into --outdir:
  reads.fq        FASTQ of reads
  truth.tsv       per-read truth (see columns below)
  answerkey.bed   founder answer-key BED (ref_chr,ref_start,ref_end,founder)

truth.tsv columns (tab-separated):
  read_id  source_genotype  transcript  strand  exon_segments  n_junctions
  error_positions  oracle_set_per_span  is_PAV
where exon_segments = "contig:refStart-refEnd;..." (0-based half-open, reference
coords) and oracle_set_per_span = "i,j|k,..." (comma-joined gamete indices per
segment, error-free, both strands), aligned with exon_segments order.
"""
import argparse
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simlib


def build_transcripts(sim, seq_of, ems, sources):
    """Return list of transcript dicts: {source, tid, strand, mrna, exons_m}
    where exons_m = [(m_start, m_end, g0, g1)] in transcription order."""
    txs = []
    for src in sources:
        em = ems[src]
        absent = set(sim["per_gt"][src]["absent_genes"])
        gseq = seq_of[src]
        for gene in sim["genes"]:
            if gene["id"] in absent:
                continue
            strand = gene["strand"]
            # canonical exons (reference coords) -> this genotype's coords
            gexons = []
            ok = True
            for (r0, r1) in gene["exons"]:
                a0 = em.ref_to_geno(r0)
                a1 = em.ref_to_geno(r1 - 1)
                if a0 is None or a1 is None:
                    ok = False
                    break
                gexons.append((a0, a1 + 1))
            if not ok:
                continue
            gexons.sort()
            order = gexons if strand == "+" else list(reversed(gexons))
            mrna_parts, exons_m = [], []
            m = 0
            for (g0, g1) in order:
                piece = gseq[g0:g1]
                if strand == "-":
                    piece = simlib.revcomp(piece)
                mrna_parts.append(piece)
                exons_m.append((m, m + len(piece), g0, g1))
                m += len(piece)
            txs.append({"source": src, "tid": "%s:%s" % (src, gene["id"]),
                        "strand": strand, "mrna": "".join(mrna_parts),
                        "exons_m": exons_m})
    return txs


def read_segments(tx, a, b):
    """For read covering mRNA [a,b), return list of segments
    (m_lo, m_hi, gen_lo, gen_hi) where [gen_lo,gen_hi) is the genotype interval
    (half-open, ascending) of that segment."""
    segs = []
    for (m_start, m_end, g0, g1) in tx["exons_m"]:
        lo = max(a, m_start)
        hi = min(b, m_end)
        if lo >= hi:
            continue
        o0, o1 = lo - m_start, hi - m_start
        if tx["strand"] == "+":
            gen_lo, gen_hi = g0 + o0, g0 + o1
        else:
            gen_lo, gen_hi = g1 - o1, g1 - o0
        segs.append((lo, hi, gen_lo, gen_hi))
    segs.sort()
    return segs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--simdir", required=True, help="sim_genomes.py output dir")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--nreads", type=int, default=20000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--error", type=float, default=0.01,
                    help="per-base substitution error rate")
    ap.add_argument("--source", default="B73",
                    help="source genotype, or 'all' for every genotype")
    ap.add_argument("--max-junctions", type=int, default=0,
                    help="reject reads spanning more than this many junctions "
                         "(0 = within-exon only, the Stage-0 default)")
    ap.add_argument("--library", choices=["sense", "antisense", "unstranded"],
                    default="sense",
                    help="read strand vs the mRNA: sense (default), antisense, or "
                         "unstranded (per-read coin flip)")
    # adversarial RNA negatives (labeled in the read_class truth column; default off)
    ap.add_argument("--intron-retention", type=float, default=0.0,
                    help="fraction of reads drawn genomic (intron kept, contiguous) "
                         "straddling a small intron")
    ap.add_argument("--chimera", type=float, default=0.0,
                    help="fraction of reads that fuse two transcript halves "
                         "(template switch; non-colinear)")
    ap.add_argument("--expr-sigma", type=float, default=1.0,
                    help="log-normal expression sigma (0 = uniform)")
    ap.add_argument("--max-tries", type=int, default=20)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    sim = json.load(open(os.path.join(a.simdir, "sim.json")))
    seqs = simlib.read_fasta(os.path.join(a.simdir, "pangenome.fa"))
    contig = sim["contig"]
    gindex = sim["gamete_index"]
    genotypes = sim["genotypes"]
    seq_of = {gt: seqs["%s_%s" % (gt, contig)] for gt in genotypes}
    ems = {gt: simlib.EditMap([tuple(b) for b in sim["per_gt"][gt]["blocks"]])
           for gt in genotypes}
    geno_seqs = {gt: seq_of[gt] for gt in genotypes}
    ref_len = sim["ref_len"]
    projections = {gt: simlib.build_ref_projection(
        [tuple(b) for b in sim["per_gt"][gt]["blocks"]], seq_of[gt], ref_len)
        for gt in genotypes}

    sources = genotypes if a.source == "all" else a.source.split(",")
    for s in sources:
        if s not in genotypes:
            sys.exit("unknown --source %s" % s)

    txs = build_transcripts(sim, seq_of, ems, sources)
    L = a.read_len
    usable = [t for t in txs if len(t["mrna"]) >= L]
    if not usable:
        sys.exit("no transcript >= read length %d" % L)

    # expression weight: log-normal per transcript, scaled by samplable positions.
    weights = []
    for t in usable:
        expr = rng.lognormvariate(0.0, a.expr_sigma) if a.expr_sigma > 0 else 1.0
        weights.append(expr * max(1, len(t["mrna"]) - L + 1))

    os.makedirs(a.outdir, exist_ok=True)
    fq = open(os.path.join(a.outdir, "reads.fq"), "w")
    tv = open(os.path.join(a.outdir, "truth.tsv"), "w")
    tv.write("read_id\tsource_genotype\ttranscript\tstrand\texon_segments\t"
             "n_junctions\terror_positions\toracle_set_per_span\tis_PAV\tread_class\n")

    sub = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    # per-transcript genomic exons (sorted) + introns, for intron-retention reads
    for t in usable:
        ge = sorted((g0, g1) for (_, _, g0, g1) in t["exons_m"])
        t["introns"] = [(ge[k][1], ge[k + 1][0]) for k in range(len(ge) - 1)]

    def span_truth(source, gen_lo, gen_hi, clean_sub):
        """(seg_str, oracle_str, is_pav) for a genotype-genomic span [gen_lo,gen_hi)."""
        ref = ems[source].geno_span_to_ref(gen_lo, gen_hi)
        if ref is None:
            oset = simlib.oracle_set(clean_sub, geno_seqs, gindex)
            return "%s:-1--1" % contig, ",".join(map(str, oset)), True
        oset = simlib.oracle_from_projection(projections, source, ref[0], ref[1], gindex)
        return "%s:%d-%d" % (contig, ref[0], ref[1]), ",".join(map(str, oset)), False

    def payload(clean, segs_with_src, njunc, source, tid, strand, rclass):
        ss, os_, pav = [], [], False
        for (src, gen_lo, gen_hi, sub_seq) in segs_with_src:
            s1, o1, p = span_truth(src, gen_lo, gen_hi, sub_seq)
            ss.append(s1); os_.append(o1); pav = pav or p
        return clean, ss, os_, njunc, pav, source, tid, strand, rclass

    def gen_normal():
        t = rng.choices(usable, weights=weights, k=1)[0]
        start = rng.randint(0, len(t["mrna"]) - L)
        segs = read_segments(t, start, start + L)
        if len(segs) - 1 > a.max_junctions:
            return None
        clean = t["mrna"][start:start + L]
        sw = [(t["source"], gl, gh, clean[ml - start:mh - start]) for (ml, mh, gl, gh) in segs]
        return payload(clean, sw, len(segs) - 1, t["source"], t["tid"], t["strand"], "normal")

    def gen_intron_retention():
        # a genomic (unspliced) read straddling a small intron -- must map contiguously,
        # not be fabricated into a junction. Only introns short enough to fit in a read.
        cand = [t for t in usable if any(0 < (b - c) <= L - 40 for (c, b) in t["introns"])]
        if not cand:
            return None
        t = rng.choice(cand)
        c, b = rng.choice([(c, b) for (c, b) in t["introns"] if 0 < (b - c) <= L - 40])
        start_g = c - rng.randint(20, L - (b - c) - 20)      # exon bases before the intron
        gseq = seq_of[t["source"]]
        if start_g < 0 or start_g + L > len(gseq):
            return None
        clean = gseq[start_g:start_g + L]                     # intron kept -> genomic-contiguous
        return payload(clean, [(t["source"], start_g, start_g + L, clean)], 0,
                       t["source"], t["tid"] + ":IR", t["strand"], "intron_retention")

    def gen_chimera():
        # two transcript halves fused (template switch): the halves are non-colinear
        # and must NOT yield a confident single-locus placement.
        if len(usable) < 2:
            return None
        tA, tB = rng.sample(usable, 2)
        h = L // 2
        aS, bS = rng.randint(0, len(tA["mrna"]) - h), rng.randint(0, len(tB["mrna"]) - (L - h))
        clean = tA["mrna"][aS:aS + h] + tB["mrna"][bS:bS + (L - h)]
        sw = []
        for (t, s, e, off) in ((tA, aS, aS + h, 0), (tB, bS, bS + (L - h), h)):
            for (ml, mh, gl, gh) in read_segments(t, s, e):
                sw.append((t["source"], gl, gh, clean[off + (ml - s):off + (mh - s)]))
        return payload(clean, sw, len(sw) - 1, tA["source"],
                       tA["tid"] + "+" + tB["tid"], ".", "chimera")

    n_made = n_skip = 0
    for i in range(a.nreads):
        u = rng.random()
        gen = (gen_chimera if u < a.chimera else
               gen_intron_retention if u < a.chimera + a.intron_retention else gen_normal)
        p = None
        for _ in range(a.max_tries):
            p = gen()
            if p is not None:
                break
        if p is None:
            n_skip += 1
            continue
        clean, seg_strs, oracle_strs, njunc, is_pav, source, tid, strand, rclass = p

        # library strand: sense = mRNA orientation; antisense = its reverse complement;
        # unstranded = a coin flip. Oracle/segments are strand-independent.
        antisense = (a.library == "antisense" or
                     (a.library == "unstranded" and rng.random() < 0.5))
        base = simlib.revcomp(clean) if antisense else clean
        r = list(base)
        errs = []
        for j, ch in enumerate(r):
            if ch in sub and rng.random() < a.error:
                r[j] = rng.choice(sub[ch])
                errs.append(j)
        read = "".join(r)

        rid = "r%06d" % i
        fq.write("@%s\n%s\n+\n%s\n" % (rid, read, "I" * len(read)))
        tv.write("\t".join([
            rid, source, tid, strand, ";".join(seg_strs), str(njunc),
            ",".join(map(str, errs)), "|".join(oracle_strs),
            "1" if is_pav else "0", rclass,
        ]) + "\n")
        n_made += 1

    fq.close()
    tv.close()

    # founder answer-key BED (Tier B). Single-source sample -> one founder over
    # the whole reference; multi-source is a placeholder (Tier B mosaics later).
    with open(os.path.join(a.outdir, "answerkey.bed"), "w") as out:
        if len(sources) == 1:
            out.write("%s\t0\t%d\t%s\n" % (contig, sim["ref_len"], sources[0]))
        else:
            out.write("%s\t0\t%d\t%s\n" % (contig, sim["ref_len"], sources[0]))

    sys.stderr.write("wrote %d reads (%d skipped) to %s\n"
                     % (n_made, n_skip, a.outdir))
    sys.stderr.write("sources=%s error=%.3f max_junctions=%d\n"
                     % (",".join(sources), a.error, a.max_junctions))


if __name__ == "__main__":
    main()
