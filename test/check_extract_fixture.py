#!/usr/bin/env python3
"""Check `ropebwt3 extract` against the simulated-pangenome fixture (simulate_pangenome.py).

  check_extract_fixture.py FIXTURE_DIR DIRECT_FA PROJ_TSV PROJ_FA [--min-proj-exact F] [--min-proj-within F N]

DIRECT_FA  `extract -f FIXTURE/direct_queries.tsv` output: every record must equal
           FIXTURE/direct_expected.fa exactly (hard failure otherwise).
PROJ_TSV   `extract --lift ... --project-only -f FIXTURE/ref_queries.tsv` output.
PROJ_FA    `extract --lift ... -f FIXTURE/ref_queries.tsv` output (FASTA, same projection).

Reports, for the reference queries: projection success rate, exact agreement of the
projected interval with the true one (true_contig/true_start/true_end), start/end
offset distributions, the same split by spans_large_indel, the sequence identity of
the projected extraction vs ref_expected.fa (global edit distance), and checks that
each PROJ_FA record is byte-identical to the fixture's own founder FASTA at the
projected coordinates (i.e. extraction itself is exact even where projection is off).
Exit status 1 if a hard check or an optional --min-* threshold fails.
"""
import argparse
import sys
from collections import Counter

import numpy as np

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fa(fn):
    d, k, buf = {}, None, []
    hdr = {}
    for line in open(fn):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if k is not None:
                d[k] = "".join(buf)
            f = line[1:].split()
            k, buf = f[0], []
            hdr[k] = f[1:]
        else:
            buf.append(line)
    if k is not None:
        d[k] = "".join(buf)
    return d, hdr


def read_tsv(fn):
    rows = []
    with open(fn) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    return rows


def edit_distance(a, b):
    """Global (Levenshtein) edit distance, Myers/Hyyro bit-parallel on Python ints."""
    if len(a) < len(b):
        a, b = b, a
    m = len(b)
    if m == 0:
        return len(a)
    peq = {}
    for i, c in enumerate(b):
        peq[c] = peq.get(c, 0) | (1 << i)
    mask = (1 << m) - 1
    hi = 1 << (m - 1)
    pv, mv, score = mask, 0, m
    for c in a:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & mask)
        mh = pv & xh
        if ph & hi:
            score += 1
        elif mh & hi:
            score -= 1
        ph = ((ph << 1) | 1) & mask
        mh = (mh << 1) & mask
        pv = mh | (~(xv | ph) & mask)
        mv = ph & xv
    return score


def dist_summary(x):
    x = np.asarray(x)
    if len(x) == 0:
        return "n=0"
    ax = np.abs(x)
    qs = np.percentile(ax, [50, 90, 99])
    return (f"n={len(x)} exact={np.mean(x == 0)*100:.1f}% |d|<=10:{np.mean(ax <= 10)*100:.1f}% "
            f"<=100:{np.mean(ax <= 100)*100:.1f}% <=1000:{np.mean(ax <= 1000)*100:.1f}% "
            f"median|d|={qs[0]:.0f} p90={qs[1]:.0f} p99={qs[2]:.0f} max={ax.max()} mean(d)={x.mean():+.1f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fixture")
    ap.add_argument("direct_fa")
    ap.add_argument("proj_tsv")
    ap.add_argument("proj_fa")
    ap.add_argument("--min-proj-exact", type=float, default=None, help="min fraction of exact intervals")
    ap.add_argument("--min-proj-within", type=float, nargs=2, default=None, metavar=("FRAC", "BP"),
                    help="min fraction with both |start|,|end| offsets <= BP")
    ap.add_argument("--no-identity", action="store_true", help="skip the edit-distance identity (slow-ish)")
    a = ap.parse_args()
    ok = True

    # ---- direct queries: exact ----
    exp, _ = read_fa(f"{a.fixture}/direct_expected.fa")
    got, _ = read_fa(a.direct_fa)
    n_eq = sum(got.get(k) == v for k, v in exp.items())
    print(f"direct: {n_eq}/{len(exp)} records identical to direct_expected.fa ({len(got)} records written)")
    if n_eq != len(exp) or len(got) != len(exp):
        bad = [k for k, v in exp.items() if got.get(k) != v][:5]
        print(f"FAIL direct: mismatching/missing e.g. {bad}")
        ok = False

    # ---- reference queries: projection ----
    q = {r["qid"]: r for r in read_tsv(f"{a.fixture}/ref_queries.tsv")}
    pr = {r["qid"]: r for r in read_tsv(a.proj_tsv)}
    rexp, _ = read_fa(f"{a.fixture}/ref_expected.fa")
    pfa, phdr = read_fa(a.proj_fa)
    asm_seq = {}
    for asm in sorted({r["assembly"] for r in q.values()}):
        asm_seq.update(read_fa(f"{a.fixture}/{asm}.fa")[0])
    n = len(q)
    n_ok = sum(1 for k in q if pr.get(k, {}).get("status", "").startswith("ok"))
    status = Counter(pr[k]["status"].split(":")[0] if k in pr else "missing" for k in q)
    print(f"projection: {n_ok}/{n} projected ({dict(status)})")
    meth = Counter(pr[k]["method"] for k in q if k in pr)
    print(f"projection methods (start,end): {dict(meth.most_common())}")
    groups = {"all": [], "spans_large_indel=0": [], "spans_large_indel=1": []}
    for k, r in q.items():
        p = pr.get(k)
        if p is None or not p["status"].startswith("ok"):
            continue
        rec = dict(qid=k, contig_ok=p["contig"] == r["true_contig"], strand=p["strand"],
                   ds=int(p["start"]) - int(r["true_start"]), de=int(p["end"]) - int(r["true_end"]),
                   method=p["method"])
        groups["all"].append(rec)
        groups[f"spans_large_indel={r['spans_large_indel']}"].append(rec)
    for g, recs in groups.items():
        nq = sum(1 for k, r in q.items() if g == "all" or g == f"spans_large_indel={r['spans_large_indel']}")
        if nq == 0:
            continue
        ex = sum(1 for x in recs if x["contig_ok"] and x["strand"] == "+" and x["ds"] == 0 and x["de"] == 0)
        print(f"[{g}] queries={nq} projected={len(recs)} contig_ok={sum(x['contig_ok'] for x in recs)} "
              f"exact_interval={ex} ({ex/nq*100:.1f}% of queries)")
        print(f"   start offset: {dist_summary([x['ds'] for x in recs if x['contig_ok']])}")
        print(f"   end offset:   {dist_summary([x['de'] for x in recs if x['contig_ok']])}")
    allr = [x for x in groups["all"] if x["contig_ok"]]
    if a.min_proj_exact is not None:
        frac = sum(1 for x in allr if x["ds"] == 0 and x["de"] == 0) / n
        if frac < a.min_proj_exact:
            print(f"FAIL: exact-interval fraction {frac:.3f} < {a.min_proj_exact}")
            ok = False
    if a.min_proj_within is not None:
        f0, bp = a.min_proj_within
        frac = sum(1 for x in allr if abs(x["ds"]) <= bp and abs(x["de"]) <= bp) / n
        if frac < f0:
            print(f"FAIL: within-{bp:.0f}bp fraction {frac:.3f} < {f0}")
            ok = False

    # ---- projected FASTA: extraction exact at the projected coordinates; identity vs truth ----
    n_fa_ok, ids, ids_by = 0, [], {"0": [], "1": []}
    for k, r in q.items():
        p = pr.get(k)
        if p is None or not p["status"].startswith("ok"):
            if k in pfa:
                print(f"FAIL: {k} has a FASTA record but no successful projection")
                ok = False
            continue
        s = pfa.get(k)
        want = asm_seq[p["contig"]][int(p["start"]):int(p["end"])]
        if p["strand"] == "-":
            want = revcomp(want)
        if s == want:
            n_fa_ok += 1
        if not a.no_identity and s is not None:
            d = edit_distance(s, rexp[k])
            ident = 1 - d / max(len(s), len(rexp[k]))
            ids.append(ident)
            ids_by[r["spans_large_indel"]].append(ident)
    print(f"projected FASTA: {n_fa_ok}/{n_ok} records equal the founder FASTA at the projected interval")
    if n_fa_ok != n_ok:
        ok = False
        print("FAIL: projected extraction differs from the founder sequence at the projected interval")
    if ids:
        for lab, v in (("all", ids), ("spans_large_indel=0", ids_by["0"]), ("spans_large_indel=1", ids_by["1"])):
            if v:
                v = np.asarray(v)
                print(f"identity vs ref_expected.fa [{lab}]: n={len(v)} mean={v.mean()*100:.3f}% "
                      f"median={np.median(v)*100:.3f}% min={v.min()*100:.2f}% ==100%: {np.mean(v == 1)*100:.1f}% "
                      f">=99%: {np.mean(v >= .99)*100:.1f}%")
    print("CHECK_EXTRACT_FIXTURE_OK" if ok else "CHECK_EXTRACT_FIXTURE_FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
