#!/usr/bin/env python3
"""Inject substitution errors into reads (models sequencing error). Read names are
preserved so the existing truth.tsv still joins. Deterministic per --seed."""
import argparse, random, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rate", type=float, default=0.01, help="per-base substitution rate")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rnd = random.Random(a.seed)
    sub = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    n = nsub = nb = 0
    with open(a.inp) as f, open(a.out, "w") as out:
        name = None
        for line in f:
            if line[0] == ">":
                name = line.rstrip()
                out.write(line)
                n += 1
            else:
                s = list(line.strip())
                for i, c in enumerate(s):
                    nb += 1
                    if c in sub and rnd.random() < a.rate:
                        s[i] = rnd.choice(sub[c]); nsub += 1
                out.write("".join(s) + "\n")
    sys.stderr.write(f"{n} reads, {nsub}/{nb} bases substituted ({100.0*nsub/nb:.2f}%)\n")


if __name__ == "__main__":
    main()
