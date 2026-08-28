#!/usr/bin/env python3
"""Shared helpers for the RNAseq PS4G simulator (Module 0). Pure stdlib.

Provides:
  - FASTA read/write
  - GFF3 write
  - Block-based genotype<->reference exact coordinate map (EditMap)
  - gamete indexing matching ropebwt3 ps4g.c: the gamete is the sequence name up
    to the first '_', gametes are numbered 0..N-1 in *sorted* unique-name order,
    so simulator indices equal the aligner's and GRITS's by construction.
  - oracle_set: exact-substring membership of a sequence across all genotypes
    (both strands) -- the ground-truth "maximal gamete set".

Coordinate model. A genotype is derived from the reference by a list of blocks:
  ('M', g0, g1, r0, r1) : genotype[g0:g1] is aligned to reference[r0:r1]
                          (g1-g0 == r1-r0; substitutions allowed, no indels
                          inside a block).
  ('I', g0, g1, r, r)   : genotype[g0:g1] is inserted, no reference image; r is
                          the reference insertion locus.
Deletions are implicit: a gap between one block's r1 and the next block's r0.
The reference genotype itself has a single identity block ('M',0,L,0,L).
"""
import gzip

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


# ---------------------------------------------------------------- FASTA / GFF

def read_fasta(path):
    """Return dict {name: uppercase_seq}. Handles .gz. Name = first token."""
    op = gzip.open if str(path).endswith(".gz") else open
    out, name, parts = {}, None, []
    with op(path, "rt") as f:
        for line in f:
            if line[0] == ">":
                if name is not None:
                    out[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        out[name] = "".join(parts).upper()
    return out


def write_fasta(path, records, width=60):
    """records: iterable of (name, seq)."""
    with open(path, "w") as out:
        for name, seq in records:
            out.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def write_gff3(path, seqid_of, genes_per_gt):
    """genes_per_gt: dict gt -> list of gene dicts with genotype coords:
       {id, strand, exons:[(g0,g1),...]}  (0-based half-open -> GFF 1-based)."""
    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for gt, genes in genes_per_gt.items():
            sid = seqid_of(gt)
            for gene in genes:
                if not gene["exons"]:
                    continue
                gs = min(g0 for g0, _ in gene["exons"])
                ge = max(g1 for _, g1 in gene["exons"])
                gid = "%s:%s" % (gt, gene["id"])
                out.write("\t".join([sid, "sim", "gene", str(gs + 1), str(ge),
                                     ".", gene["strand"], ".",
                                     "ID=%s" % gid]) + "\n")
                out.write("\t".join([sid, "sim", "mRNA", str(gs + 1), str(ge),
                                     ".", gene["strand"], ".",
                                     "ID=%s.t1;Parent=%s" % (gid, gid)]) + "\n")
                for k, (g0, g1) in enumerate(gene["exons"]):
                    out.write("\t".join([sid, "sim", "exon", str(g0 + 1), str(g1),
                                         ".", gene["strand"], ".",
                                         "ID=%s.t1.e%d;Parent=%s.t1" % (gid, k, gid)]) + "\n")


# ------------------------------------------------------------------ gametes

def gamete_of(seqname):
    """Sample name = sequence name up to the first '_'  (B73_chr1 -> B73)."""
    i = seqname.find("_")
    return seqname if i < 0 else seqname[:i]


def contig_of(seqname):
    """refContig = the part after the first '_'  (B73_chr1 -> chr1)."""
    i = seqname.find("_")
    return seqname if i < 0 else seqname[i + 1:]


def gamete_index(names):
    """names: iterable of gamete (sample) names. Returns {name: index} numbered
    0..N-1 in sorted unique-name order -- identical to ropebwt3 ps4g.c."""
    uniq = sorted(set(names))
    return {n: i for i, n in enumerate(uniq)}


# ------------------------------------------------------------------ EditMap

class EditMap:
    """Exact genotype<->reference coordinate map for one genotype."""

    def __init__(self, blocks):
        # blocks: list of (kind, g0, g1, r0, r1); sorted by g0.
        self.blocks = sorted(blocks, key=lambda b: b[1])
        self._g0 = [b[1] for b in self.blocks]

    def geno_to_ref(self, gpos):
        """Reference base aligned to genotype position gpos, or None (inserted)."""
        import bisect
        j = bisect.bisect_right(self._g0, gpos) - 1
        if j < 0:
            return None
        kind, g0, g1, r0, r1 = self.blocks[j]
        if gpos >= g1:
            return None
        if kind == "M":
            return r0 + (gpos - g0)
        return None  # inserted base: no reference image

    def geno_span_to_ref(self, ga, gb):
        """Map genotype span [ga,gb) to a reference interval (r_lo, r_hi) using
        only the aligned ('M') bases it covers, or None if it covers none
        (fully inserted / PAV). r_hi is exclusive."""
        rs = []
        for kind, g0, g1, r0, r1 in self.blocks:
            if kind != "M":
                continue
            lo = max(ga, g0)
            hi = min(gb, g1)
            if lo < hi:
                rs.append(r0 + (lo - g0))
                rs.append(r0 + (hi - g0) - 1)
        if not rs:
            return None
        return min(rs), max(rs) + 1

    def ref_to_geno(self, rpos):
        """Genotype position aligned to reference base rpos, or None (deleted)."""
        for kind, g0, g1, r0, r1 in self.blocks:
            if kind == "M" and r0 <= rpos < r1:
                return g0 + (rpos - r0)
        return None


# ------------------------------------------------------------------ oracle

def oracle_set(seq, geno_seqs, gindex):
    """Set of gamete indices whose genotype sequence contains `seq` (or its
    reverse complement) as an exact substring, searched anywhere. `geno_seqs`
    maps gamete name -> full genotype sequence. O(genome) per genotype; used for
    the rare inserted/PAV segment that has no reference image."""
    rc = revcomp(seq)
    out = set()
    for name, gseq in geno_seqs.items():
        if seq in gseq or rc in gseq:
            out.add(gindex[name])
    return sorted(out)


def build_ref_projection(blocks, gseq, ref_len, gap="."):
    """Project a genotype onto reference coordinates: a string of length
    ref_len whose r-th char is the genotype's base aligned to reference base r,
    or `gap` where the genotype has deleted that reference base."""
    arr = [gap] * ref_len
    for (kind, g0, g1, r0, r1) in blocks:
        if kind == "M":
            arr[r0:r1] = gseq[g0:g1]
    return "".join(arr)


def oracle_from_projection(projections, src, r_lo, r_hi, gindex):
    """Oracle set at an *orthologous* reference interval: the gamete indices
    whose reference projection is identical to the source's over [r_lo, r_hi).
    `projections` maps gamete name -> build_ref_projection() string. For a
    colinear exon read this equals the exact-substring oracle (a 150 bp match
    elsewhere in the pangenome is astronomically unlikely), while being O(L) per
    genotype instead of O(genome)."""
    key = projections[src][r_lo:r_hi]
    out = []
    for name, proj in projections.items():
        if proj[r_lo:r_hi] == key:
            out.append(gindex[name])
    return sorted(out)
