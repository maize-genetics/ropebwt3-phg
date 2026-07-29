#!/usr/bin/env python3
"""E4 prototype, step 1: sample anchor k-mers along the carrier genomes.

Each anchor is a short window taken at a fixed stride from a carrier chromosome.
Running `ropebwt3 mem -p` on these anchors yields, for every anchor that is
shared with the reference, BOTH its carrier coordinate and its B73 coordinate
(same locus) -> the (carrier_pos <-> B73_pos) pairs used to build a per-chromosome
liftover. The anchor name encodes its true carrier origin: CARRIER|chr|pos|len.

Usage:
    gen_anchors.py --len 100 --stride 10000 --out anchors.fa \
        B97:B97.chr.fa.gz Ki3:Ki3.chr.fa.gz CML247:CML247.chr.fa.gz
"""
import argparse, gzip, sys


def read_chrs(path):
    chroms = {}
    name = None
    parts = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line[0] == ">":
                if name is not None:
                    chroms[name] = "".join(parts).upper()
                full = line[1:].split()[0]
                name = full.split("_", 1)[1] if "_" in full else full
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        chroms[name] = "".join(parts).upper()
    return chroms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--len", type=int, default=100)
    ap.add_argument("--stride", type=int, default=10000)
    ap.add_argument("--out", required=True)
    ap.add_argument("inputs", nargs="+", help="CARRIER:path.fa.gz")
    a = ap.parse_args()
    L = a.len
    n = 0
    with open(a.out, "w") as out:
        for spec in a.inputs:
            carrier, path = spec.split(":", 1)
            chroms = read_chrs(path)
            for c, seq in chroms.items():
                for pos in range(0, len(seq) - L, a.stride):
                    w = seq[pos:pos + L]
                    if "N" in w:
                        continue
                    out.write(f">{carrier}|{c}|{pos}|{L}\n{w}\n")
                    n += 1
    sys.stderr.write(f"wrote {n} anchors to {a.out}\n")


if __name__ == "__main__":
    main()
