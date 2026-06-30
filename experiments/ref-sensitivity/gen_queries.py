#!/usr/bin/env python3
"""Draw random fixed-length queries from renamed (chr-only) NAM genomes.

Each input is a gzipped FASTA whose records are named "<LINE>_chr1".."<LINE>_chr10".
For every line we sample N windows of length L, uniformly over the concatenated
chromosome length (so longer chromosomes get proportionally more queries), rejecting
any window that contains an 'N'. The query name encodes ground truth:

    >{LINE}|{chr}|{pos0}|{len}        e.g.  B97|chr5|10234567|150

Also writes a TSV truth table: qname<TAB>line<TAB>chr<TAB>pos0<TAB>len.

Usage:
    gen_queries.py --n 25000 --len 150 --seed 7 \
        --out queries.fa --truth truth.tsv \
        B73:data/B73.chr.fa.gz B97:data/B97.chr.fa.gz ...
"""
import argparse, gzip, random, sys


def read_chroms(path):
    """Return dict {chr_suffix: sequence} for one renamed genome (chr-only)."""
    chroms, name, parts = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == ">":
                if name is not None:
                    chroms[name] = "".join(parts).upper()
                # header like "B97_chr5" -> keep "chr5"
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
    ap.add_argument("--n", type=int, required=True, help="queries per genome")
    ap.add_argument("--len", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--max-tries", type=int, default=50)
    ap.add_argument("inputs", nargs="+", help="LINE:path.gz pairs")
    a = ap.parse_args()

    rng = random.Random(a.seed)
    L = a.len
    fa = open(a.out, "w")
    tv = open(a.truth, "w")
    tv.write("qname\tline\tchr\tpos0\tlen\n")

    for spec in a.inputs:
        line, path = spec.split(":", 1)
        sys.stderr.write(f"[{line}] reading {path}\n"); sys.stderr.flush()
        chroms = read_chroms(path)
        names = sorted(chroms)
        # weight chromosome choice by (length - L), i.e. number of valid start positions
        weights = [max(0, len(chroms[c]) - L) for c in names]
        total = sum(weights)
        sys.stderr.write(f"[{line}] {len(names)} chroms, {total/1e6:.0f} Mbp samplable\n")
        made = 0
        while made < a.n:
            c = rng.choices(names, weights=weights, k=1)[0]
            seq = chroms[c]
            hi = len(seq) - L
            ok = False
            for _ in range(a.max_tries):
                p = rng.randint(0, hi)
                w = seq[p:p + L]
                if "N" not in w:
                    ok = True
                    break
            if not ok:
                continue
            qname = f"{line}|{c}|{p}|{L}"
            fa.write(f">{qname}\n{w}\n")
            tv.write(f"{qname}\t{line}\t{c}\t{p}\t{L}\n")
            made += 1
        sys.stderr.write(f"[{line}] wrote {made}\n")

    fa.close(); tv.close()
    sys.stderr.write("done\n")


if __name__ == "__main__":
    main()
