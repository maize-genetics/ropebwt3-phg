#!/usr/bin/env python3
"""Boundaries placed inside small deletions / exactly at insertions (simulated-pangenome fixture).

  check_extract_boundaries.py FIXTURE_DIR RB3 IDX.fmd LIFT [--flank 3000] [--min-del F] [--min-ins F]
                              [-- extra extract args]

From FIXTURE/events.tsv: every deletion of >= 2 bp gives a boundary x strictly inside it, every
insertion a boundary x exactly at it. Each x is split into [x-flank, x) and [x, x+flank) in the same
founder. Truth (maps/*.npy): true end of the left piece = homolog of the last present base before x,
+1; true start of the right piece = homolog of the first present base from x on (equal for a
deletion; the insertion in between for an insertion).
Deletions: "exact" if both ends are at the truth, "equivalent" if both pieces meet at a junction of an
equivalent placement of the same deletion (shifted by s with identical reference bases), which the
sequence cannot distinguish. Insertions: exact if both ends are at the truth.
Exit 1 if (exact+equivalent)/n for deletions < --min-del, or exact/n for insertions < --min-ins.
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


def read_fa(fn):
    d, k = {}, None
    for line in open(fn):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]
            d[k] = []
        else:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fixture")
    ap.add_argument("rb3")
    ap.add_argument("fmd")
    ap.add_argument("lift")
    ap.add_argument("--flank", type=int, default=3000)
    ap.add_argument("--min-del", type=float, default=None)
    ap.add_argument("--min-ins", type=float, default=None)
    ap.add_argument("extra", nargs="*")
    a = ap.parse_args()
    fx = Path(a.fixture)
    ref = read_fa(fx / "B73.fa")
    lens = {k: len(v) for k, v in ref.items()}
    cases = []
    for line in list(open(fx / "events.tsv"))[1:]:
        asm, chrom, rp, typ, rl, al = line.rstrip("\n").split("\t")
        rp, rl = int(rp), int(rl)
        rc = f"B73_{chrom}"
        if typ == "del" and rl >= 2:
            x = rp + rl // 2 if rl > 2 else rp + 1
        elif typ == "ins":
            x = rp
        else:
            continue
        if x - a.flank < 0 or x + a.flank > lens[rc]:
            continue
        cases.append((typ, asm, chrom, rc, x))
    with tempfile.TemporaryDirectory() as td:
        qf = Path(td) / "q.tsv"
        with open(qf, "w") as f:
            f.write("qid\tref_contig\tref_start\tref_end\tassembly\n")
            for k, (typ, asm, chrom, rc, x) in enumerate(cases):
                f.write(f"b{k}L\t{rc}\t{x - a.flank}\t{x}\t{asm}\nb{k}R\t{rc}\t{x}\t{x + a.flank}\t{asm}\n")
        out = subprocess.run([a.rb3, "extract", f"--lift={a.lift}", "--ref-prefix=B73", "--project-only"] + a.extra
                             + ["-f", str(qf), a.fmd], capture_output=True, text=True, check=True).stdout
    hdr, *lines = out.rstrip("\n").split("\n")
    hdr = hdr.split("\t")
    pr = {d["qid"]: d for d in (dict(zip(hdr, l.split("\t"))) for l in lines)}
    maps, n = {}, {"del": 0, "ins": 0}
    ex, eq, meet, fail = dict(n), dict(n), dict(n), dict(n)
    for k, (typ, asm, chrom, rc, x) in enumerate(cases):
        n[typ] += 1
        L, R = pr[f"b{k}L"], pr[f"b{k}R"]
        if not (L["status"].startswith("ok") and R["status"].startswith("ok")) or L["contig"] != R["contig"]:
            fail[typ] += 1
            continue
        key = f"{asm}_{chrom}"
        if key not in maps:
            maps[key] = np.load(fx / "maps" / f"{key}.npy")
        mp = maps[key]
        pl = x - 1
        while pl >= 0 and mp[pl] < 0:
            pl -= 1
        pf = x
        while mp[pf] < 0:
            pf += 1
        te, ts = int(mp[pl]) + 1, int(mp[pf])
        le, rs = int(L["end"]), int(R["start"])
        meet[typ] += (rs - le) == (ts - te)
        if le == te and rs == ts:
            ex[typ] += 1
        elif typ == "del" and le == rs:
            s, d0, d1, r = le - te, pl + 1, pf, ref[rc]
            eq[typ] += (s > 0 and r[d0:d0 + s] == r[d1:d1 + s]) or (s < 0 and r[d0 + s:d0] == r[d1 + s:d1])
    for t in ("del", "ins"):
        print(f"{t}: n={n[t]} failed={fail[t]} pieces meet as in the truth={meet[t]} exact={ex[t]}"
              + (f" equivalent={eq[t]}" if t == "del" else ""))
    ok = True
    if a.min_del is not None and n["del"] and (ex["del"] + eq["del"]) / n["del"] < a.min_del:
        ok = False
    if a.min_ins is not None and n["ins"] and ex["ins"] / n["ins"] < a.min_ins:
        ok = False
    print("CHECK_EXTRACT_BOUNDARIES_OK" if ok else "CHECK_EXTRACT_BOUNDARIES_FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
