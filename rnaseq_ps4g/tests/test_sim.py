#!/usr/bin/env python3
"""Tests for the Module-0 simulator: parameterized divergence targets, exact
coordinate-map integrity, exon conservation (SVs confined to intron/intergenic),
and the recall==1.0 oracle self-consistency invariant.

Run:  pytest rnaseq_ps4g/tests/
"""
import itertools
import json
import os
import subprocess
import sys

try:
    import pytest
except ImportError:                       # allow running as a plain script
    pytest = None

HERE = os.path.dirname(os.path.abspath(__file__))
SIMDIR = os.path.join(HERE, "..", "sim")
sys.path.insert(0, SIMDIR)
import simlib  # noqa: E402


def _run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def build_sim(base):
    """Build one small synthetic pangenome + a zero-error read set."""
    gdir = os.path.join(base, "genomes")
    rdir = os.path.join(base, "reads")
    _run([sys.executable, os.path.join(SIMDIR, "sim_genomes.py"),
          "--outdir", gdir, "--seed", "7", "--ngenes", "5",
          "--intergenic-min", "8000", "--intergenic-max", "15000"])
    _run([sys.executable, os.path.join(SIMDIR, "sim_rnaseq.py"),
          "--simdir", gdir, "--outdir", rdir,
          "--nreads", "3000", "--seed", "11", "--error", "0",
          "--source", "all", "--max-junctions", "0"])
    j = json.load(open(os.path.join(gdir, "sim.json")))
    seqs = simlib.read_fasta(os.path.join(gdir, "pangenome.fa"))
    ems = {gt: simlib.EditMap([tuple(b) for b in j["per_gt"][gt]["blocks"]])
           for gt in j["genotypes"]}
    seq_of = {gt: seqs["%s_%s" % (gt, j["contig"])] for gt in j["genotypes"]}
    return {"j": j, "seq_of": seq_of, "ems": ems, "gdir": gdir, "rdir": rdir}


if pytest is not None:
    @pytest.fixture(scope="module")
    def sim(tmp_path_factory):
        return build_sim(str(tmp_path_factory.mktemp("sim")))


def _base_at_ref(sim, gt, r):
    g = sim["ems"][gt].ref_to_geno(r)
    return None if g is None else sim["seq_of"][gt][g]


def _positions_by_kind(sim, kind):
    out = []
    for (r0, r1, k) in sim["j"]["regions"]:
        if k == kind:
            out.extend(range(r0, r1))
    return out


def test_gamete_index_sorted_rule(sim):
    names = sim["j"]["genotypes"]
    assert sim["j"]["gamete_index"] == simlib.gamete_index(names)
    # equals sorted order
    assert sim["j"]["gamete_index"] == {n: i for i, n in enumerate(sorted(names))}


def test_exon_pairwise_divergence(sim):
    gts = sim["j"]["genotypes"]
    pos = _positions_by_kind(sim, "exon")
    pos = pos[::max(1, len(pos) // 5000)]
    diff = den = 0
    for a, b in itertools.combinations(gts, 2):
        for r in pos:
            ba, bb = _base_at_ref(sim, a, r), _base_at_ref(sim, b, r)
            if ba is not None and bb is not None:
                den += 1
                diff += (ba != bb)
    div = diff / den
    assert 0.006 <= div <= 0.014, "exon pairwise SNP divergence %.4f !~ 0.01" % div


def test_intron_intergenic_shared_half(sim):
    gts = sim["j"]["genotypes"]
    pos = _positions_by_kind(sim, "intron") + _positions_by_kind(sim, "intergenic")
    pos = pos[::max(1, len(pos) // 8000)]
    shared = total = 0
    for a, b in itertools.combinations(gts, 2):
        for r in pos:
            total += 1
            ba, bb = _base_at_ref(sim, a, r), _base_at_ref(sim, b, r)
            if ba is not None and bb is not None and ba == bb:
                shared += 1
    frac = shared / total
    assert 0.40 <= frac <= 0.60, "intron+intergenic shared %.4f !~ 0.5" % frac


def test_editmap_roundtrip(sim):
    for gt in sim["j"]["genotypes"]:
        em = sim["ems"][gt]
        for (kind, g0, g1, r0, r1) in em.blocks:
            if kind != "M":
                continue
            for g in (g0, (g0 + g1) // 2, g1 - 1):
                r = em.geno_to_ref(g)
                assert r is not None and em.ref_to_geno(r) == g


def test_exons_conserved_not_partially_disrupted(sim):
    """Each canonical exon is either fully present or fully absent (PAV) in a
    genotype -- structural variation is confined to introns/intergenic."""
    for gt in sim["j"]["genotypes"]:
        em = sim["ems"][gt]
        for gene in sim["j"]["genes"]:
            for (r0, r1) in gene["exons"]:
                present = [em.ref_to_geno(r) is not None for r in range(r0, r1)]
                assert all(present) or not any(present), \
                    "exon %d-%d partially disrupted in %s" % (r0, r1, gt)


def test_oracle_self_consistency(sim):
    """The recall==1.0 invariant: for an error-free read, the source genotype is
    in the oracle set of every segment. Tool-independent ground-truth check."""
    gidx = sim["j"]["gamete_index"]
    n = 0
    with open(os.path.join(sim["rdir"], "truth.tsv")) as f:
        f.readline()
        for line in f:
            F = line.rstrip("\n").split("\t")
            src, errpos, oracle_str = F[1], F[6], F[7]
            assert errpos == "", "error-0 run produced a read with errors"
            for seg in oracle_str.split("|"):
                idxs = {int(x) for x in seg.split(",") if x != ""}
                assert gidx[src] in idxs, "source %s missing from oracle %s" % (src, seg)
            n += 1
    assert n > 0


def _main():
    """Run the suite without pytest (repo plain-script convention)."""
    import tempfile
    base = tempfile.mkdtemp(prefix="test_sim_")
    s = build_sim(base)
    tests = [test_gamete_index_sorted_rule, test_exon_pairwise_divergence,
             test_intron_intergenic_shared_half, test_editmap_roundtrip,
             test_exons_conserved_not_partially_disrupted,
             test_oracle_self_consistency]
    fail = 0
    for t in tests:
        try:
            t(s)
            print("PASS", t.__name__)
        except AssertionError as e:
            fail += 1
            print("FAIL", t.__name__, "--", e)
    print("%d/%d passed" % (len(tests) - fail, len(tests)))
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    _main()
