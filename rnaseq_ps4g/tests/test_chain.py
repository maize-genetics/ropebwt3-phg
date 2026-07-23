#!/usr/bin/env python3
"""Unit tests for the colinear chaining + strict intersection core
(chain/chain_prototype.py: chain_and_intersect), independent of the index.

Covers: strict in-order chaining, negative-strand genes, junction segmentation,
whole-read set intersection (spurious removal + junction linkage), rejection of
non-colinear (chimeric) and out-of-order SMEMs.

Run:  python3 rnaseq_ps4g/tests/test_chain.py   (or via pytest)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "chain"))
import chain_prototype as cp


def smem(qs, qe, refpos, strand, gametes, count=None, refs=None):
    r = refs if refs is not None else [(refpos, strand)]
    return {"qs": qs, "qe": qe, "count": count if count is not None else len(gametes),
            "gametes": set(gametes), "refs": r, "ref": r[0] if r else None,
            "informative": True}


def test_within_exon_plus_intersects():
    # two colinear '+' SMEMs, different sets -> one segment, intersection
    segs, g, flag = cp.chain_and_intersect(
        [smem(0, 60, 1000, "+", {0, 1, 2}), smem(63, 150, 1063, "+", {0, 4})], 30)
    assert flag == "ok" and g == {0}, (flag, g)
    assert segs == [(1000, 1150)], segs


def test_within_exon_minus():
    # negative strand: reference DECREASES as query increases
    segs, g, flag = cp.chain_and_intersect(
        [smem(0, 29, 215800, "-", {0, 1, 2}), smem(30, 150, 215679, "-", {0, 4})], 30)
    assert flag == "ok" and g == {0}, (flag, g)
    assert segs and segs[0][0] == 215679, segs


def test_junction_minus_two_segments_and_linkage():
    # spliced '-' read: exon1 (all 5) + exon2 (missing gamete 3) at an intron-gap.
    exon1 = smem(0, 21, 181733, "-", {0, 1, 2, 3, 4})
    exon2 = smem(118, 150, 181493, "-", {0, 1, 2, 4})
    segs, g, flag = cp.chain_and_intersect([exon1, exon2], 30)
    assert flag == "ok" and g == {0, 1, 2, 4}, (flag, g)
    assert len(segs) == 2, segs                       # two exon segments (intron gap)


def test_large_intron_placed_not_rejected():
    # exons 5 kb apart (large intron) must still chain (the gap penalty allows it
    # because it is the only 2-anchor option) -> placed as two segments.
    exon1 = smem(0, 70, 100000, "+", {0, 1, 2, 4})
    exon2 = smem(70, 150, 105000, "+", {0, 1, 2, 4})
    segs, g, flag = cp.chain_and_intersect([exon1, exon2], 30, max_intron=200000)
    assert flag == "ok" and g == {0, 1, 2, 4}, (flag, g)
    assert len(segs) == 2, segs


def test_chimera_penalized_prefers_compact_chain():
    # a long chimeric SMEM at a distant locus loses to a compact 2-exon chain of the
    # SAME anchor count because its unexplained jump is penalized.
    exon1 = smem(0, 21, 1000, "+", {0, 1, 2, 3, 4})
    exon2 = smem(118, 150, 1200, "+", {0, 1, 2, 4})
    chimera = smem(20, 117, 40000, "+", {3})          # far -> heavy gap penalty
    segs, g, flag = cp.chain_and_intersect([exon1, exon2, chimera], 30)
    assert flag == "ok" and g == {0, 1, 2, 4}, (flag, g)


def test_gtag_splice_check():
    # canonical GT..AG intron between two exons -> chained (2 segments); the same
    # gap over non-canonical reference -> the link is rejected (only one exon).
    ref = "A" * 1000 + "GT" + "C" * 196 + "AG" + "A" * 1000   # intron [1000,1200)
    ej = smem(0, 70, 930, "+", {0, 1})      # ref [930,1000) -> donor at 1000
    ei = smem(70, 150, 1200, "+", {0, 1})   # ref [1200,1280) -> acceptor before 1200
    segs, g, flag = cp.chain_and_intersect([ej, ei], 30, max_intron=2000, refseq=ref)
    assert flag == "ok" and len(segs) == 2, (flag, segs)         # canonical -> spliced
    segs, g, flag = cp.chain_and_intersect([ej, ei], 30, max_intron=2000, refseq="A" * 3000)
    assert len(segs) == 1, segs                                  # non-canonical -> not chained


def test_tandem_dup_ambiguous():
    # a SMEM with two reference occurrences (tandem duplication) -> the read maps to
    # two loci -> flagged ambiguous, not placed at an arbitrary one.
    dup = smem(0, 150, 1000, "+", {0, 1}, refs=[(1000, "+"), (1715, "+")])
    segs, g, flag = cp.chain_and_intersect([dup], 30)
    assert flag == "ambiguous", (flag, segs, g)


def test_out_of_order_dropped():
    a = smem(0, 50, 1000, "+", {0, 1})
    bad = smem(60, 110, 500, "+", {2})                # ref < previous -> out of order
    c = smem(120, 150, 1100, "+", {0, 3})
    segs, g, flag = cp.chain_and_intersect([a, bad, c], 30)
    assert flag == "ok" and g == {0}, (flag, g)       # bad SMEM excluded


def _main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    fail = 0
    for t in tests:
        try:
            t()
            print("PASS", t.__name__)
        except AssertionError as e:
            fail += 1
            print("FAIL", t.__name__, "--", e)
    print("%d/%d passed" % (len(tests) - fail, len(tests)))
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    _main()
