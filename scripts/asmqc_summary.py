#!/usr/bin/env python3
"""One-row-per-line flag summary of `ropebwt3 asmqc` / `asmqc-fa` output.

Usage:
  asmqc_summary.py --qc asmqc.tsv [--chains chains.tsv] [--anchorwave phg_metrics.csv]
                   [--out-tsv summary.tsv] [--out-json summary.json]

Flags, per reference chromosome (numbered chromosomes only), two tiers:
  check : frac_mapped < --check-mapped, or inverted / other-chromosome / scaffold bp >= --check-mb
  note  : frac_mapped < --note-mapped,  or those bp >= --note-mb
Inversions are batch-aware: a large inverted chain whose breakpoints match (within --shared-tol-mb) an
inverted chain in another line of the same batch is treated as a shared, likely real, inversion
polymorphism (reference carries the other arrangement) and reported as a 'note' naming those lines,
not a 'check' (needs --chains). Screen lines together so shared inversions can be recognised.
plus a genome-wide 'check' when the line's whole-genome mapped fraction < --genome-mapped, and a
'note' for divergence estimates above the batch median + 3 MAD. Defaults were set from the 24 NAM
founders (background p90: 4.9 Mb inverted, 2.3 Mb other-chromosome, 3.8 Mb scaffold).

--chains adds where flagged sequence went (largest chains per flag); --anchorwave adds the PHGv2
AnchorWave metrics of the same line when present, for side-by-side review.
Standard library only.
"""
import argparse, csv, json, re, statistics, sys


def chrom_num(name):
    m = re.search(r'chr0*(\d+)$', name, re.I)
    return int(m.group(1)) if m else None


def load_qc(fn):
    lines = {}
    for r in csv.DictReader(open(fn), delimiter='\t'):
        lines.setdefault(r['taxon'], {})[r['chrom']] = r
    return lines


def load_chains(fn):
    by = {}
    if not fn:
        return by
    for r in csv.DictReader(open(fn), delimiter='\t'):
        taxon = r['asm_seq'].split('_', 1)[0]
        by.setdefault((taxon, r['ref_seq'], r['class']), []).append(r)
    return by


def load_anchorwave(fn):
    aw = {}
    if not fn:
        return aw
    rows = list(csv.reader(open(fn)))
    hdr = rows[0]
    for r in rows[1:]:
        if len(r) < 5:
            continue
        # PHG metrics layout: name, method, chrom|ALL, '', metrics...; B73Ref rows put chrom in column 2
        chrom = r[2] if r[2] else r[1]
        if not (chrom == 'ALL' or chrom.isdigit()):
            continue
        d = {h: r[i] for i, h in enumerate(hdr) if i >= 4 and h}
        aw.setdefault(r[0], {})[chrom] = d
    return aw


def mb(x):
    return round(float(x) / 1e6, 2)


def base_metrics(r):
    """Base-level (asmqc-fa) columns, AnchorWave names; empty when the run had no base level."""
    if 'numSNPs' not in r or (r.get('percentRefAligned') in (None, '', '0.000000')):
        return {}
    return {'snps': int(r['numSNPs']), 'indels': int(r['numIndels']), 'aligned': round(float(r['percentRefAligned']), 4),
            'deleted': round(float(r['percentRefDeleted']), 4), 'bp_ins_mb': mb(r['numBasesInserted']),
            'bp_del_mb': mb(r['numBasesDeleted']), 'median_indel': int(r['medianIndelSize']), 'largest_indel_mb': mb(r['largestIndel'])}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--qc', required=True)
    ap.add_argument('--chains')
    ap.add_argument('--anchorwave')
    ap.add_argument('--out-tsv')
    ap.add_argument('--out-json')
    ap.add_argument('--check-mapped', type=float, default=0.90)
    ap.add_argument('--note-mapped', type=float, default=0.95)
    ap.add_argument('--check-mb', type=float, default=10.0)
    ap.add_argument('--note-mb', type=float, default=5.0)
    ap.add_argument('--genome-mapped', type=float, default=0.97)
    ap.add_argument('--shared-tol-mb', type=float, default=1.0)
    a = ap.parse_args()

    qc, chains, aw = load_qc(a.qc), load_chains(a.chains), load_anchorwave(a.anchorwave)
    divs = [float(v['ALL']['div_est']) for v in qc.values() if 'ALL' in v]
    med = statistics.median(divs) if divs else 0.0
    mad = statistics.median([abs(d - med) for d in divs]) if divs else 0.0

    # large inverted chains of every line, for the shared-breakpoint test
    tol = a.shared_tol_mb * 1e6
    inv_all = [(t, c, int(r['ref_start']), int(r['ref_end'])) for (t, c, cls), lst in chains.items() if cls == 'syn_inv'
               for r in lst if int(r['ref_end']) - int(r['ref_start']) >= a.note_mb * 1e6]

    def shared_with(taxon, c, st, en):
        return sorted({t for t, cc, s2, e2 in inv_all if t != taxon and cc == c and abs(s2 - st) <= tol and abs(e2 - en) <= tol})

    out = []
    for taxon in sorted(qc):
        rows = qc[taxon]
        g = rows.get('ALL')
        chroms = sorted((c for c in rows if c != 'ALL' and chrom_num(c) is not None), key=chrom_num)
        flags, per = [], []
        for c in chroms:
            r = rows[c]
            n = chrom_num(c)
            fm = float(r['frac_mapped'])
            entry = {'chrom': n, 'ref_len_mb': mb(r['ref_len']), 'mapped': round(fm, 4),
                     'inv_mb': mb(r['inv_bp']), 'other_mb': mb(r['other_chr_bp']), 'scaf_mb': mb(r['scaf_bp']),
                     'ins_mb': mb(r['ins_bp_lb']), 'del_mb': mb(r['del_bp_lb']), 'div_est': float(r['div_est'])}
            entry.update(base_metrics(r))
            awc = aw.get(taxon, {}).get(str(n))
            if awc:
                entry['aw_deleted'] = round(float(awc['percentRefDeleted']), 4)
                entry['aw_indels'] = int(float(awc['numIndels']))
                entry['aw_mapped'] = round(float(awc['percentRefMapped']), 4)
                entry['aw_aligned'] = round(float(awc['percentRefAligned']), 4)
                entry['aw_snps'] = int(float(awc['numSNPs']))
                entry['aw_largest_del_mb'] = mb(awc['largestDelSize'])
            per.append(entry)
            if r['n_pts'] == '0':
                flags.append({'sev': 'check', 'chrom': n, 'kind': 'no anchors', 'value': 0,
                              'text': f'chr{n}: no anchors at all (chromosome missing or misnamed?)'})
                continue
            if fm < a.check_mapped or fm < a.note_mapped:
                sev = 'check' if fm < a.check_mapped else 'note'
                flags.append({'sev': sev, 'chrom': n, 'kind': 'low mapped', 'value': round(fm, 4),
                              'text': f'chr{n}: only {fm*100:.1f}% of B73 in syntenic chains'})
            for key, cls, label in (('inv_bp', 'syn_inv', 'inverted'), ('other_chr_bp', 'other_chr', 'on other chromosomes'),
                                    ('scaf_bp', 'scaffold', 'on unplaced scaffolds')):
                v = float(r[key]) / 1e6
                if v >= a.note_mb:
                    sev = 'check' if v >= a.check_mb else 'note'
                    top = sorted(chains.get((taxon, c, cls), []), key=lambda x: -(int(x['ref_end']) - int(x['ref_start'])))[:3]
                    where = [{'ref': f"chr{n}:{mb(t['ref_start'])}-{mb(t['ref_end'])} Mb",
                              'asm': f"{t['asm_seq']}:{mb(t['asm_start'])}-{mb(t['asm_end'])} Mb", 'strand': t['strand'],
                              'span_mb': mb(int(t['ref_end']) - int(t['ref_start']))} for t in top]
                    text = f'chr{n}: {v:.1f} Mb {label}'
                    if cls == 'syn_inv' and top:
                        # bp of this line's large inverted chains that recur in other lines at the same breakpoints
                        sh_bp, sh_lines = 0, set()
                        for t in top:
                            o = shared_with(taxon, c, int(t['ref_start']), int(t['ref_end']))
                            if o:
                                sh_bp += int(t['ref_end']) - int(t['ref_start'])
                                sh_lines.update(o)
                                next(w for w in where if w['asm'].startswith(t['asm_seq'] + ':') and w['ref'].endswith(f"{mb(t['ref_end'])} Mb"))['shared_with'] = o
                        if sh_bp >= 0.5 * float(r[key]):
                            sev = 'note'
                            text += f" at breakpoints shared with {', '.join(sorted(sh_lines))}: likely a real inversion polymorphism"
                    flags.append({'sev': sev, 'chrom': n, 'kind': label, 'value': round(v, 2), 'text': text, 'where': where})
        gm = float(g['frac_mapped']) if g else 0.0
        if g and gm < a.genome_mapped:
            flags.insert(0, {'sev': 'check', 'chrom': None, 'kind': 'genome mapped', 'value': round(gm, 4),
                             'text': f'genome: only {gm*100:.1f}% of B73 in syntenic chains'})
        dv = float(g['div_est']) if g else 0.0
        if g and mad > 0 and dv > med + 3 * mad:
            flags.append({'sev': 'note', 'chrom': None, 'kind': 'divergent', 'value': dv,
                          'text': f'divergence estimate {dv*1e3:.2f}/kb, above batch median {med*1e3:.2f} + 3 MAD'})
        n_check = sum(f['sev'] == 'check' for f in flags)
        n_note = sum(f['sev'] == 'note' for f in flags)
        rec = {'line': taxon, 'status': 'check' if n_check else ('note' if n_note else 'clean'),
               'n_check': n_check, 'n_note': n_note,
               'genome': {'mapped': round(gm, 4), 'inv_mb': mb(g['inv_bp']) if g else 0, 'other_mb': mb(g['other_chr_bp']) if g else 0,
                          'scaf_mb': mb(g['scaf_bp']) if g else 0, 'ins_mb': mb(g['ins_bp_lb']) if g else 0,
                          'del_mb': mb(g['del_bp_lb']) if g else 0, 'div_est': dv, **(base_metrics(g) if g else {})},
               'flags': flags, 'chroms': per}
        awg = aw.get(taxon, {}).get('ALL')
        if awg:
            rec['aw_genome'] = {'mapped': round(float(awg['percentRefMapped']), 4), 'aligned': round(float(awg['percentRefAligned']), 4),
                                'deleted': round(float(awg['percentRefDeleted']), 4), 'median_indel': int(float(awg['medianIndelSize'])),
                                'snps': int(float(awg['numSNPs'])), 'indels': int(float(awg['numIndels'])),
                                'ins_mb': mb(awg['numBasesInserted']), 'del_mb': mb(awg['numBasesDeleted'])}
        out.append(rec)
    order = {'check': 0, 'note': 1, 'clean': 2}
    out.sort(key=lambda r: (order[r['status']], -r['n_check'], -r['n_note'], r['line']))

    if a.out_json:
        json.dump({'thresholds': vars(a), 'div_median': med, 'div_mad': mad, 'lines': out}, open(a.out_json, 'w'), indent=1)
    fo = open(a.out_tsv, 'w') if a.out_tsv else sys.stdout
    fo.write('line\tstatus\tn_check\tn_note\tgenome_mapped\tinv_mb\tother_chr_mb\tscaf_mb\tins_mb_lb\tdel_mb_lb\tdiv_est_per_kb\tSNPs\taligned\tbp_ins_mb\tbp_del_mb\tflags\n')
    for r in out:
        g = r['genome']
        fo.write(f"{r['line']}\t{r['status']}\t{r['n_check']}\t{r['n_note']}\t{g['mapped']:.4f}\t{g['inv_mb']}\t{g['other_mb']}\t{g['scaf_mb']}"
                 f"\t{g['ins_mb']}\t{g['del_mb']}\t{g['div_est']*1e3:.2f}\t{g.get('snps','')}\t{g.get('aligned','')}\t{g.get('bp_ins_mb','')}\t{g.get('bp_del_mb','')}"
                 f"\t{'; '.join(f['sev'].upper()+': '+f['text'] for f in r['flags'])}\n")


if __name__ == '__main__':
    main()
