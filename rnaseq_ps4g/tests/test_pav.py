#!/usr/bin/env python3
"""End-to-end tests for PAV breakpoint anchoring (`ropebwt3 chain --lift`).

Simulates a pangenome containing structural insertions that are ABSENT from the
reference and present in a known subset of assemblies, then checks that reads
lying inside one are placed at its breakpoint and attributed to exactly the
assemblies that carry it.

This needs the simulator's --pav-insert-* options: the churn's own insertions are
--churn-max (70bp) at most, below the read length, so without them no simulated
read is ever assembly-only and this path is unreachable.

Requires a built ../../ropebwt3; skipped if absent.

Run:  pytest rnaseq_ps4g/tests/test_pav.py
"""
import collections
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
RB = os.path.abspath(os.path.join(HERE, "..", "..", "ropebwt3"))

READ_LEN = 150
STEP = 75
GRID = 500          # --pav-grid: reported-position resolution (independent of --pav-agree)
TILED_GT = "CML247"  # the assembly reads are tiled from


def _run(cmd, **kw):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, **kw)


def build_pav_sim(base):
    """Simulate a pangenome with structural insertions, index it, build the
    liftover, and tile reads across one carrying assembly."""
    # Gentle churn on purpose. At the default rates the intergenic regions retain
    # only ~half their bases, so most reads have no exact reference match and ~97%
    # of rows come out as pav:, burying the real insertions in noise.
    _run([sys.executable, os.path.join(SIMDIR, "sim_genomes.py"),
          "--outdir", base, "--seed", "7", "--ngenes", "5",
          "--intergenic-min", "8000", "--intergenic-max", "15000",
          "--pav-insert-n", "3", "--pav-insert-bp", "2000", "--pav-insert-share", "0.6",
          "--p-del", "0.04", "--p-ins", "0.01", "--p-dup", "0.01",
          "--intergenic-snp", "0.004", "--intron-snp", "0.004"])

    idx = os.path.join(base, "idx.fmd")
    _run([RB, "build", "-d", "-o", idx, os.path.join(base, "pangenome.fa")])
    _run([RB, "ssa", "-s8", "-o", idx + ".ssa", idx])
    # sequence names/lengths sidecar
    import gzip
    lens = []
    for line in open(os.path.join(base, "pangenome.fa")):
        if line.startswith(">"):
            lens.append([line[1:].split()[0], 0])
        else:
            lens[-1][1] += len(line.strip())
    with gzip.open(idx + ".len.gz", "wt") as fh:
        for name, ln in lens:
            fh.write("%s\t%d\n" % (name, ln))
    # small -k/-s: the defaults (100/2000) leave too few anchors at this genome size,
    # and breakpoint resolution is bounded by the anchor stride.
    _run([RB, "lift", "--ref-prefix=B73", "-k", "50", "-s", "200",
          "-o", os.path.join(base, "idx.lift"), idx, os.path.join(base, "B73.fa")])

    seq = "".join(l.strip() for l in open(os.path.join(base, "%s.fa" % TILED_GT))
                  if not l.startswith(">"))
    reads = os.path.join(base, "tiles.fa")
    with open(reads, "w") as fh:
        for i in range(0, len(seq) - READ_LEN, STEP):
            fh.write(">t|%d\n%s\n" % (i, seq[i:i + READ_LEN]))
    return {"dir": base, "idx": idx, "lift": os.path.join(base, "idx.lift"),
            "reads": reads, "sim": json.load(open(os.path.join(base, "sim.json")))}


def run_chain(env, grid=GRID, extra=()):
    out = os.path.join(env["dir"], "out%s.tsv" % grid)
    cmd = [RB, "chain", "--ref-prefix=B73", "--lift=" + env["lift"],
           "--pav-grid=%d" % grid, "-l", "31", "-t", "4",
           *extra, env["idx"], env["reads"]]
    with open(out, "w") as fh:
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.DEVNULL)
    rows = [l.split("\t") for l in open(out).read().splitlines()[1:]]
    return [{"read": r[0], "contig": r[1], "pos": int(r[2]),
             "set": frozenset(int(x) for x in r[3].split(",")),
             "cls": int(r[4]) if len(r) > 4 else None} for r in rows if len(r) >= 4]


def _truth(env):
    gi = {g: i for i, g in enumerate(env["sim"]["genotypes"])}
    return gi, env["sim"]["pav_inserts"]


# --------------------------------------------------------------------------- #

def test_insertions_recovered_with_exact_assembly_set(env):
    """Every insertion the tiled assembly carries is placed near its true
    breakpoint, attributed to exactly the assemblies that carry it."""
    gi, inserts = _truth(env)
    rows = [r for r in run_chain(env) if r["contig"].startswith("pav:")]
    assert rows, "no pav rows emitted at all"
    mine = [d for d in inserts if TILED_GT in d["carriers"]]
    assert mine, "simulation produced no insertion carried by " + TILED_GT
    for d in mine:
        want = frozenset(gi[c] for c in d["carriers"])
        near = [r for r in rows
                if abs(r["pos"] - d["locus"]) <= 2 * GRID + 500 and r["set"] == want]
        assert near, ("insertion at %d (carriers %s) not recovered with set %s"
                      % (d["locus"], d["carriers"], sorted(want)))
        assert len(near) >= 5, ("only %d reads support the insertion at %d"
                                % (len(near), d["locus"]))


def test_non_carried_insertion_is_absent(env):
    """An insertion the tiled assembly does NOT carry must not be attributed to
    it -- the negative control that makes the positive result meaningful."""
    gi, inserts = _truth(env)
    rows = run_chain(env)
    others = [d for d in inserts if TILED_GT not in d["carriers"]]
    if not others:
        return                              # nothing to check in this draw
    for d in others:
        want = frozenset(gi[c] for c in d["carriers"])
        hits = [r for r in rows if r["set"] == want
                and abs(r["pos"] - d["locus"]) <= 2 * GRID + 500]
        assert not hits, ("reads from %s were attributed to the %s-only insertion at %d"
                          % (TILED_GT, ",".join(d["carriers"]), d["locus"]))


def test_pav_rows_never_contain_the_reference(env):
    """A pav: row means the sequence is absent from the reference, so assembly
    index 0 (B73) can never appear in one."""
    for r in run_chain(env):
        if r["contig"].startswith("pav:"):
            assert 0 not in r["set"], "reference assembly present in a pav row: %r" % (r,)


def test_insertion_rows_are_classified_as_insertions(env):
    """Reads inside sequence the reference lacks are class 1 (insertion), not
    class 0 (diverged)."""
    gi, inserts = _truth(env)
    rows = run_chain(env)
    mine = [d for d in inserts if TILED_GT in d["carriers"]]
    for d in mine:
        want = frozenset(gi[c] for c in d["carriers"])
        near = [r for r in rows
                if abs(r["pos"] - d["locus"]) <= 2 * GRID + 500 and r["set"] == want]
        assert near and all(r["cls"] == 1 for r in near), \
            "insertion at %d not classified as an insertion" % d["locus"]


def test_breakpoint_is_consistent_across_reads(env):
    """All reads inside one insertion should agree on a single breakpoint once
    snapped to the grid -- otherwise one insertion's evidence fragments across
    bins and each bin is individually weaker."""
    gi, inserts = _truth(env)
    rows = run_chain(env)
    for d in [x for x in inserts if TILED_GT in x["carriers"]]:
        want = frozenset(gi[c] for c in d["carriers"])
        near = [r for r in rows
                if abs(r["pos"] - d["locus"]) <= 2 * GRID + 500 and r["set"] == want]
        assert near
        spread = max(r["pos"] for r in near) - min(r["pos"] for r in near)
        assert spread <= 2 * GRID, \
            ("insertion at %d scatters over %d bp (grid %d)" % (d["locus"], spread, GRID))


def test_min_len_floor_suppresses_pav_rows(env):
    """--pav-min-len above the read length must suppress every pav row while
    leaving reference-anchored rows alone."""
    rows = run_chain(env, extra=["--pav-min-len=200"])
    assert not [r for r in rows if r["contig"].startswith("pav:")]
    assert [r for r in rows if not r["contig"].startswith("pav:")], \
        "reference-anchored rows were suppressed too"


def test_grid_changes_resolution_not_membership(env):
    """--pav-grid sets the resolution of the reported position and nothing else.
    It must not change WHICH reads are emitted, their assembly sets, or their
    class -- only where they are reported. (It used to also set the agreement
    threshold, which made --pav-grid=0 the strictest setting rather than the
    loosest; the two are now separate knobs.)"""
    snapped = run_chain(env, grid=GRID)
    exact = run_chain(env, grid=0)
    key = lambda rows: [(r["read"], r["set"], r["cls"]) for r in rows]
    assert key(snapped) == key(exact), "grid changed which reads/sets were emitted"
    assert [r["pos"] for r in snapped] != [r["pos"] for r in exact], \
        "grid did not change the reported positions at all"


def test_agree_suppresses_when_tightened(env):
    """--pav-agree is the suppression knob: an assembly's projections of the seed
    must fall within it. Tightening it can only remove rows, never add them."""
    loose = [r for r in run_chain(env, extra=["--pav-agree=50000"])
             if r["contig"].startswith("pav:")]
    tight = [r for r in run_chain(env, extra=["--pav-agree=30"])
             if r["contig"].startswith("pav:")]
    assert len(tight) <= len(loose)


# --------------------------------------------------------------------------- #

if pytest is not None:
    @pytest.fixture(scope="module")
    def env(tmp_path_factory):
        if not os.access(RB, os.X_OK):
            pytest.skip("ropebwt3 not built at %s; run make" % RB)
        return build_pav_sim(str(tmp_path_factory.mktemp("pav")))
else:
    def _main():
        import tempfile
        if not os.access(RB, os.X_OK):
            print("SKIP: ropebwt3 not built at %s" % RB); return 0
        with tempfile.TemporaryDirectory() as td:
            e = build_pav_sim(td)
            fails = 0
            for name, fn in sorted(globals().items()):
                if name.startswith("test_") and callable(fn):
                    try:
                        fn(e); print("ok   %s" % name)
                    except AssertionError as ex:
                        fails += 1; print("FAIL %s: %s" % (name, ex))
            print("%d failed" % fails)
            return 1 if fails else 0

    if __name__ == "__main__":
        sys.exit(_main())
