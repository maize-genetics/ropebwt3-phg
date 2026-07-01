#!/usr/bin/env python3
"""Tile each read into fixed-length k-mers. Each k-mer is emitted as its own FASTA
record named READNAME#OFFSET so the placer can regroup k-mers by their read and
know each k-mer's position within it."""
import argparse, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=75)
    ap.add_argument("--step", type=int, default=25)
    a = ap.parse_args()
    nreads = nk = 0
    with open(a.inp) as f, open(a.out, "w") as out:
        name = None; seq = []
        def flush():
            nonlocal nk
            if name is None:
                return
            s = "".join(seq)
            last = -1
            for off in range(0, max(1, len(s) - a.k + 1), a.step):
                last = off
                out.write(f">{name}#{off}\n{s[off:off+a.k]}\n"); nk += 1
            if len(s) >= a.k and last != len(s) - a.k:      # ensure the 3' end is covered
                off = len(s) - a.k
                out.write(f">{name}#{off}\n{s[off:off+a.k]}\n"); nk += 1
        for line in f:
            if line[0] == ">":
                flush()
                name = line[1:].split()[0]; seq = []; nreads += 1
            else:
                seq.append(line.strip())
        flush()
    sys.stderr.write(f"{nreads} reads -> {nk} k-mers (k={a.k}, step={a.step})\n")


if __name__ == "__main__":
    main()
