#!/usr/bin/env python3
"""Boundary consistency of `ropebwt3 extract --lift` on the simulated-pangenome fixture.

  check_extract_junctions.py FIXTURE_DIR RB3 IDX.fmd LIFT [--min-consistent F] [-- extra extract args]

Every reference query [s, e) of FIXTURE/ref_queries.tsv is split at m = (s+e)//2 and both halves are
projected into the same founder (an X -> X "switch"). The stitched junction offset is
right.start - left.end. The truth (maps/*.npy) gives the true end of [s, m) and the true start of
[m, e); their difference is the founder sequence that belongs to neither piece (bases inserted
exactly at m; 0 for no indel at m and for a deletion containing m). The projection is consistent
when its offset equals that truth gap. Exit 1 if the consistent fraction is below --min-consistent.
"""
import argparse
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fixture")
    ap.add_argument("rb3")
    ap.add_argument("fmd")
    ap.add_argument("lift")
    ap.add_argument("--min-consistent", type=float, default=None)
    ap.add_argument("extra", nargs="*")
    a = ap.parse_args()
    fx = Path(a.fixture)
    rows = [l.rstrip("\n").split("\t") for l in open(fx / "ref_queries.tsv")][1:]
    with tempfile.TemporaryDirectory() as td:
        qf = Path(td) / "q.tsv"
        with open(qf, "w") as f:
            f.write("qid\tref_contig\tref_start\tref_end\tassembly\n")
            for r in rows:
                s, e = int(r[2]), int(r[3])
                m = (s + e) // 2
                f.write(f"{r[0]}|L\t{r[1]}\t{s}\t{m}\t{r[4]}\n{r[0]}|R\t{r[1]}\t{m}\t{e}\t{r[4]}\n")
        out = subprocess.run([a.rb3, "extract", f"--lift={a.lift}", "--ref-prefix=B73", "--project-only"] + a.extra
                             + ["-f", str(qf), a.fmd], capture_output=True, text=True, check=True).stdout
    hdr, *lines = out.rstrip("\n").split("\n")
    hdr = hdr.split("\t")
    pr = {}
    for l in lines:
        d = dict(zip(hdr, l.split("\t")))
        pr[d["qid"]] = d
    n = n_cons = n_zero = n_fail = 0
    bad = Counter()
    examples = []
    for r in rows:
        s, e = int(r[2]), int(r[3])
        m = (s + e) // 2
        L, R = pr.get(r[0] + "|L"), pr.get(r[0] + "|R")
        if L is None or R is None or not L["status"].startswith("ok") or not R["status"].startswith("ok") \
                or L["contig"] != R["contig"] or L["strand"] != R["strand"]:
            n_fail += 1
            continue
        mp = np.load(fx / "maps" / f"{r[4]}_{r[1].split('_', 1)[1]}.npy")
        pl = m - 1
        while pl >= 0 and mp[pl] < 0:
            pl -= 1
        pf = m
        while pf < len(mp) - 1 and mp[pf] < 0:
            pf += 1
        truth_gap = int(mp[pf]) - (int(mp[pl]) + 1)
        off = int(R["start"]) - int(L["end"]) if L["strand"] == "+" else int(L["start"]) - int(R["end"])
        n += 1
        n_zero += off == 0
        if off == truth_gap:
            n_cons += 1
        else:
            d = off - truth_gap
            bad["<=10" if abs(d) <= 10 else "<=100" if abs(d) <= 100 else ">100"] += 1
            if len(examples) < 10:
                examples.append((r[0], r[4], m, off, truth_gap, L.get("snap", ""), R.get("snap", "")))
    frac = n_cons / n if n else 0
    print(f"junctions: {n} scored ({n_fail} with a failed/split projection); consistent with the truth gap: "
          f"{n_cons} ({frac*100:.1f}%); offset exactly 0: {n_zero}; inconsistent by |d| {dict(bad)}")
    for x in examples:
        print("  inconsistent:", x)
    if a.min_consistent is not None and frac < a.min_consistent:
        print("CHECK_EXTRACT_JUNCTIONS_FAILED")
        sys.exit(1)
    print("CHECK_EXTRACT_JUNCTIONS_OK")


if __name__ == "__main__":
    main()
