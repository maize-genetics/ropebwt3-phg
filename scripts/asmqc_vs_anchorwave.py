#!/usr/bin/env python3
"""Agreement between `ropebwt3 asmqc-fa` output and PHGv2 AnchorWave assembly metrics.

Usage:
  asmqc_vs_anchorwave.py --qc asmqc.tsv --anchorwave phg_metrics.csv [--json out.json]

For each metric: Pearson r across lines (genome totals over chr1-10), Pearson r across
line x chromosome rows, and the median asmqc/AnchorWave ratio across lines. Lines are matched
by name; chromosomes by number. Standard library only.
"""
import argparse, csv, json, math, re, statistics, sys

METRICS = [  # label, asmqc column, AnchorWave column, kind ('frac' = fraction of reference)
    ('Mapped fraction', 'frac_mapped', 'percentRefMapped', 'frac'),
    ('Aligned fraction', 'percentRefAligned', 'percentRefAligned', 'frac'),
    ('Deleted fraction', 'percentRefDeleted', 'percentRefDeleted', 'frac'),
    ('Inserted bp', 'numBasesInserted', 'numBasesInserted', 'count'),
    ('Deleted bp', 'numBasesDeleted', 'numBasesDeleted', 'count'),
    ('Indels', 'numIndels', 'numIndels', 'count'),
    ('SNPs', 'numSNPs', 'numSNPs', 'count'),
]


def pearson(x, y):
    if len(x) < 3:
        return float('nan')
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sx = math.sqrt(sum((a - mx) ** 2 for a in x)); sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy) if sx and sy else float('nan')


def chrom_num(name):
    m = re.search(r'chr0*(\d+)$', name, re.I)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--qc', required=True); ap.add_argument('--anchorwave', required=True); ap.add_argument('--json')
    a = ap.parse_args()
    aw = {}
    rows = list(csv.reader(open(a.anchorwave))); hdr = rows[0]
    for r in rows[1:]:
        if len(r) > 4 and r[2].isdigit() and r[1]:
            aw[(r[0], int(r[2]))] = {h: r[i] for i, h in enumerate(hdr) if i >= 4 and h}
    qc = {}
    for r in csv.DictReader(open(a.qc), delimiter='\t'):
        n = chrom_num(r['chrom'])
        if n is not None and 1 <= n <= 10:
            qc[(r['taxon'], n)] = r
    keys = sorted(k for k in qc if k in aw)
    lines = sorted({k[0] for k in keys if all((k[0], c) in aw and (k[0], c) in qc for c in range(1, 11))})
    out = {'n_lines': len(lines), 'n_rows': len(keys), 'metrics': []}

    def value(src, k, col, kind):
        d = src[k]
        return float(d[col])

    for label, qcol, acol, kind in METRICS:
        rq, ra = [], []
        for k in keys:
            rq.append(value(qc, k, qcol, kind)); ra.append(value(aw, k, acol, kind))
        gq, ga = [], []
        for t in lines:
            ks = [(t, c) for c in range(1, 11)]
            if kind == 'frac':  # length-weighted genome fraction
                L = [float(aw[k]['totalBases']) for k in ks]
                gq.append(sum(value(qc, k, qcol, kind) * l for k, l in zip(ks, L)) / sum(L))
                ga.append(sum(value(aw, k, acol, kind) * l for k, l in zip(ks, L)) / sum(L))
            else:
                gq.append(sum(value(qc, k, qcol, kind) for k in ks)); ga.append(sum(value(aw, k, acol, kind) for k in ks))
        ratios = [x / y for x, y in zip(gq, ga) if y]
        out['metrics'].append({'label': label, 'aw': acol, 'r_lines': round(pearson(gq, ga), 3), 'r_rows': round(pearson(rq, ra), 3),
                               'ratio': f"{statistics.median(ratios):.2f}x" if ratios else '-',
                               'per_line': {t: [round(x, 6), round(y, 6)] for t, x, y in zip(lines, gq, ga)}})
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)
    print(f"{out['n_lines']} lines, {out['n_rows']} line x chromosome rows")
    print('metric\tAnchorWave\tr_lines\tr_rows\tmedian_ratio')
    for m in out['metrics']:
        print(f"{m['label']}\t{m['aw']}\t{m['r_lines']}\t{m['r_rows']}\t{m['ratio']}")


if __name__ == '__main__':
    main()
