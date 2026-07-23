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


def smem(qs, qe, refpos, strand, gametes, count=None):
    return {"qs": qs, "qe": qe, "count": count if count is not None else len(gametes),
            "gametes": set(gametes), "ref": (refpos, strand),
            "informative": True}


def test_within_exon_plus_intersects():
    # two colinear '+' SMEMs, different sets -> one segment, intersection
    segs, g = cp.chain_and_intersect(
        [smem(0, 60, 1000, "+", {0, 1, 2}), smem(63, 150, 1063, "+", {0, 4})], 30)
    assert g == {0}, g
    assert segs == [(1000, 1150)], segs


def test_within_exon_minus():
    # negative strand: reference DECREASES as query increases
    segs, g = cp.chain_and_intersect(
        [smem(0, 29, 215800, "-", {0, 1, 2}), smem(30, 150, 215679, "-", {0, 4})], 30)
    assert g == {0}, g
    # segments are (min_ref, max_ref) forward-strand intervals
    assert segs and segs[0][0] == 215679, segs


def test_junction_minus_two_segments_and_linkage():
    # spliced '-' read: exon1 (all 5) + exon2 (missing gamete 3) at an intron-gap.
    # whole-read intersection excludes 3 (it lacks exon2); both exons emitted.
    exon1 = smem(0, 21, 181733, "-", {0, 1, 2, 3, 4})
    exon2 = smem(118, 150, 181493, "-", {0, 1, 2, 4})
    segs, g = cp.chain_and_intersect([exon1, exon2], 30)
    assert g == {0, 1, 2, 4}, g
    assert len(segs) == 2, segs                       # two exon segments (intron gap)


def test_chimeric_noncolinear_excluded():
    # a chimeric junction-spanning SMEM with a reference hit at a NON-colinear pos
    # (breaks query<->ref order) must be dropped, not intersected (else it would
    # empty the set and drop the source).
    exon1 = smem(0, 21, 1000, "+", {0, 1, 2, 3, 4})
    exon2 = smem(118, 150, 1200, "+", {0, 1, 2, 4})
    chimera = smem(20, 117, 5000, "+", {3})           # ref jumps forward then back
    segs, g = cp.chain_and_intersect([exon1, exon2, chimera], 30)
    assert 0 in g and g == {0, 1, 2, 4}, g            # source survives; chimera gone


def test_out_of_order_dropped():
    # an out-of-order SMEM (ref goes backward on '+' strand) is not chained
    a = smem(0, 50, 1000, "+", {0, 1})
    bad = smem(60, 110, 500, "+", {2})                # ref < previous -> out of order
    c = smem(120, 150, 1100, "+", {0, 3})
    segs, g = cp.chain_and_intersect([a, bad, c], 30)
    assert 2 not in g, g                              # bad SMEM excluded
    assert g == {0}, g


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
