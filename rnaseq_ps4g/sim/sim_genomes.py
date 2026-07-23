#!/usr/bin/env python3
"""Module 0, part 1: build a small synthetic multi-genotype pangenome + GFF.

Organizing principle (HANDOFF section 4.1): exons are the conserved anchor;
introns and intergenic space are the churn.

  - one ancestral/reference genotype (default B73), a few hundred kb, with 5-6
    genes (exon/intron structure);
  - each other genotype is derived from the reference by an exact edit script:
      * exons: substitutions only (~conserved, colinear);
      * introns + intergenic: substitutions + insertions/deletions/duplications
        so a genotype pair shares only ~half of these bases base-for-base;
      * PAV: whole genes are occasionally deleted in a non-reference genotype.

Every genotype carries an exact block map to the reference (simlib.EditMap), so
downstream truth (read reference coordinates, oracle sets, answer-key BED) is
exact by construction. Outputs into --outdir:

  <GT>.fa                per-genotype FASTA (sequence named <GT>_chr1)
  pangenome.fa           all genotypes concatenated (index/aligner input)
  annot.gff3             genes/exons per genotype, in each genotype's coords
  gametes.tsv            gameteIndex<TAB>sampleName  (sorted-name rule)
  sim.json               machine-readable truth: blocks, regions, genes, PAV
"""
import argparse
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simlib

BASES = "ACGT"


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(rng, s, rate):
    """Substitute each base to a different base with probability `rate`.
    Substituted (non-ancestral) bases are emitted in lowercase so mutations are
    visible in the FASTA; ancestral bases keep the reference's uppercase.
    (ropebwt3 folds case, and read_fasta uppercases, so this is cosmetic only.)"""
    if rate <= 0:
        return s
    sub = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    out = list(s)
    for i, c in enumerate(out):
        if c in sub and rng.random() < rate:
            out[i] = rng.choice(sub[c]).lower()
    return "".join(out)


def build_reference(rng, a):
    """Return (ref_seq, genes, regions). genes: list of dicts with reference
    coords {id, strand, exons:[(r0,r1)], introns:[(r0,r1)]}. regions: list of
    (r0, r1, kind) tiling the reference with kind in exon/intron/intergenic."""
    seq = []
    genes = []
    regions = []
    pos = 0

    def add(kind, length):
        nonlocal pos
        return add_seq(kind, rand_seq(rng, length))

    def add_seq(kind, s):
        """Append a given sequence as a region (used to emit an identical tandem copy)."""
        nonlocal pos
        seq.append(s)
        r0 = pos
        pos += len(s)
        regions.append((r0, pos, kind))
        return r0, pos

    tandem = set(range(a.tandem_dup))                 # gene indices to duplicate in tandem
    large_intron_gene = a.ngenes - 1 if a.large_intron_kb > 0 else -1

    add("intergenic", rng.randint(a.intergenic_min, a.intergenic_max))
    for gi in range(a.ngenes):
        nex = rng.randint(a.min_exons, a.max_exons)
        strand = rng.choice("+-")
        # build the gene's sequence pieces once (so a tandem copy is identical)
        pieces = []                                   # list of (kind, sequence)
        for ei in range(nex):
            # ensure the first exon of gene 0 is long (clean within-exon reads);
            # sprinkle one short internal exon to exercise 2-junction reads later.
            if ei == 0:
                elen = rng.randint(max(a.read_len, 200), 400)
            elif 0 < ei < nex - 1 and rng.random() < 0.25:
                elen = rng.randint(70, 100)
            else:
                elen = rng.randint(120, 350)
            pieces.append(("exon", rand_seq(rng, elen)))
            if ei < nex - 1:
                # one designated gene gets a large first intron (splice gap > max_ref_span)
                if gi == large_intron_gene and ei == 0:
                    ilen = int(a.large_intron_kb * 1000)
                else:
                    ilen = rng.randint(a.intron_min, a.intron_max)
                pieces.append(("intron", rand_seq(rng, ilen)))
        # emit the gene, then an identical adjacent copy if it is a tandem duplicate
        for copy in range(2 if gi in tandem else 1):
            gid = "gene%d" % gi if copy == 0 else "gene%ddup" % gi
            exons, introns = [], []
            for kind, s in pieces:
                r0, r1 = add_seq(kind, s)
                (exons if kind == "exon" else introns).append((r0, r1))
            genes.append({"id": gid, "strand": strand, "exons": exons, "introns": introns})
        add("intergenic", rng.randint(a.intergenic_min, a.intergenic_max))

    return "".join(seq), genes, regions


def derive_genotype(rng, refseq, regions, genes, gt, a):
    """Build one non-reference genotype from the reference. Returns
    (geno_seq, blocks, absent_gene_ids)."""
    # decide PAV: never drop gene0 (keep a fully-conserved anchor gene).
    absent = set()
    for g in genes[1:]:
        if rng.random() < a.pav_rate:
            absent.add(g["id"])
    exon_of_absent = set()
    intron_of_absent = set()
    for g in genes:
        if g["id"] in absent:
            for e in g["exons"]:
                exon_of_absent.add(e)
            for it in g["introns"]:
                intron_of_absent.add(it)

    out = []          # genotype sequence pieces
    blocks = []       # (kind, g0, g1, r0, r1)
    gpos = 0

    def emit_match(r0, r1, rate):
        nonlocal gpos
        s = mutate(rng, refseq[r0:r1], rate)
        out.append(s)
        blocks.append(("M", gpos, gpos + len(s), r0, r1))
        gpos += len(s)

    def emit_insert(r_locus, s):
        nonlocal gpos
        if not s:
            return
        out.append(s)
        blocks.append(("I", gpos, gpos + len(s), r_locus, r_locus))
        gpos += len(s)

    for (r0, r1, kind) in regions:
        if kind == "exon":
            if (r0, r1) in exon_of_absent:
                continue  # PAV deletion: no genotype bases, ref gap
            emit_match(r0, r1, a.exon_snp)
        elif kind == "intron" and (r0, r1) in intron_of_absent:
            continue      # intron of a deleted gene
        else:
            rate = a.intron_snp if kind == "intron" else a.intergenic_snp
            _churn_region(rng, refseq, r0, r1, rate, a, emit_match, emit_insert)

    return "".join(out), blocks, sorted(absent)


def _churn_region(rng, refseq, r0, r1, snp, a, emit_match, emit_insert):
    """Walk a reference region left->right applying deletions, insertions and
    duplications so that (with an independent draw per genotype) a genotype pair
    retains only ~half of the region base-for-base."""
    reflen = len(refseq)
    r = r0
    while r < r1:
        u = rng.random()
        if u < a.p_del:
            d = min(rng.randint(a.churn_min, a.churn_max), r1 - r)
            r += d                                   # deletion: skip ref bases
            continue
        if u < a.p_del + a.p_ins:
            k = rng.randint(a.churn_min, a.churn_max)
            emit_insert(r, rand_seq(rng, k).lower())        # novel insertion (non-ancestral)
            continue
        if u < a.p_del + a.p_ins + a.p_dup:
            k = rng.randint(a.churn_min, a.churn_max)
            src = rng.randint(0, max(0, reflen - k))
            emit_insert(r, refseq[src:src + k].lower())     # duplication, non-ancestral at this locus
            continue
        run = min(rng.randint(a.churn_min, a.churn_max), r1 - r)
        emit_match(r, r + run, snp)                  # retained run (substituted)
        r += run


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--genotypes", default="B73,CML247,Ki3,Mo17,W22")
    ap.add_argument("--ref-name", default="B73")
    ap.add_argument("--contig", default="chr1")
    ap.add_argument("--ngenes", type=int, default=6)
    ap.add_argument("--min-exons", type=int, default=3)
    ap.add_argument("--max-exons", type=int, default=6)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--intron-min", type=int, default=80)
    ap.add_argument("--intron-max", type=int, default=150)
    ap.add_argument("--intergenic-min", type=int, default=30000)
    ap.add_argument("--intergenic-max", type=int, default=55000)
    ap.add_argument("--exon-snp", type=float, default=0.007,
                    help="per-genotype exon substitution rate (pairwise ~2x)")
    ap.add_argument("--intron-snp", type=float, default=0.01)
    ap.add_argument("--intergenic-snp", type=float, default=0.01)
    ap.add_argument("--pav-rate", type=float, default=0.12,
                    help="prob a gene (not gene0) is absent in a non-ref genotype")
    # adversarial fixtures (default off; turn on to stress chaining ties / gap cost)
    ap.add_argument("--tandem-dup", type=int, default=0,
                    help="duplicate the first N genes in tandem (identical adjacent "
                         "copy -> placement ties)")
    ap.add_argument("--large-intron-kb", type=float, default=0,
                    help="give the last gene one intron this many kb (splice gap "
                         "beyond the chainer's max-ref-span)")
    # churn operation probabilities (per step); remainder is a retained run.
    ap.add_argument("--p-del", type=float, default=0.30)
    ap.add_argument("--p-ins", type=float, default=0.10)
    ap.add_argument("--p-dup", type=float, default=0.06)
    ap.add_argument("--churn-min", type=int, default=15)
    ap.add_argument("--churn-max", type=int, default=70)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    genotypes = a.genotypes.split(",")
    if a.ref_name not in genotypes:
        sys.exit("--ref-name %s must be in --genotypes" % a.ref_name)
    os.makedirs(a.outdir, exist_ok=True)

    refseq, genes, regions = build_reference(rng, a)
    reflen = len(refseq)

    seqs = {}          # gt -> genotype sequence
    per_gt = {}        # gt -> {blocks, absent_genes}
    for gt in genotypes:
        if gt == a.ref_name:
            seqs[gt] = refseq
            per_gt[gt] = {"blocks": [("M", 0, reflen, 0, reflen)], "absent_genes": []}
        else:
            gseq, blocks, absent = derive_genotype(rng, refseq, regions, genes, gt, a)
            seqs[gt] = gseq
            per_gt[gt] = {"blocks": blocks, "absent_genes": absent}

    seqid_of = lambda g: "%s_%s" % (g, a.contig)
    gindex = simlib.gamete_index(genotypes)

    # per-genotype genes in genotype coordinates (for the GFF), via EditMap.
    genes_per_gt = {}
    for gt in genotypes:
        em = simlib.EditMap(per_gt[gt]["blocks"])
        absent = set(per_gt[gt]["absent_genes"])
        glist = []
        for g in genes:
            if g["id"] in absent:
                continue
            gexons = []
            ok = True
            for (r0, r1) in g["exons"]:
                a0 = em.ref_to_geno(r0)
                a1 = em.ref_to_geno(r1 - 1)
                if a0 is None or a1 is None:
                    ok = False
                    break
                gexons.append((a0, a1 + 1))
            if ok:
                glist.append({"id": g["id"], "strand": g["strand"], "exons": gexons})
        genes_per_gt[gt] = glist

    # ---- write outputs
    for gt in genotypes:
        simlib.write_fasta(os.path.join(a.outdir, "%s.fa" % gt),
                           [(seqid_of(gt), seqs[gt])])
    simlib.write_fasta(os.path.join(a.outdir, "pangenome.fa"),
                       [(seqid_of(gt), seqs[gt]) for gt in genotypes])
    simlib.write_gff3(os.path.join(a.outdir, "annot.gff3"), seqid_of, genes_per_gt)

    with open(os.path.join(a.outdir, "gametes.tsv"), "w") as out:
        for name in sorted(gindex, key=lambda n: gindex[n]):
            out.write("%d\t%s\n" % (gindex[name], name))

    sim = {
        "ref": a.ref_name,
        "contig": a.contig,
        "genotypes": genotypes,
        "gamete_index": gindex,
        "ref_len": reflen,
        "read_len": a.read_len,
        "genes": [{"id": g["id"], "strand": g["strand"],
                   "exons": g["exons"], "introns": g["introns"]} for g in genes],
        "regions": regions,
        "per_gt": per_gt,
    }
    with open(os.path.join(a.outdir, "sim.json"), "w") as out:
        json.dump(sim, out)

    # ---- brief report
    sys.stderr.write("reference %s: %d bp, %d genes, %d regions\n"
                     % (a.ref_name, reflen, len(genes), len(regions)))
    for gt in genotypes:
        sys.stderr.write("  %-8s len=%7d  PAV-absent=%s\n"
                         % (gt, len(seqs[gt]), per_gt[gt]["absent_genes"] or "-"))
    sys.stderr.write("wrote FASTAs, pangenome.fa, annot.gff3, gametes.tsv, sim.json to %s\n"
                     % a.outdir)


if __name__ == "__main__":
    main()
