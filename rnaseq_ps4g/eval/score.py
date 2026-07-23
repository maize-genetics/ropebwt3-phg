#!/usr/bin/env python3
"""Module 0, Tier-A scoring harness (HANDOFF section 4.4 / 6).

Scores an aligner run against the simulator's exact oracle truth. Consumes:
  --truth    truth.tsv from sim_rnaseq.py
  --refmap   ropebwt3 refmap per-read table (its stdout)
  --ps4g     the PS4G file refmap wrote (--ps4g)

Scored from the individual's perspective: reads are sampled from one gamete
(homozygous) now, a heterozygous individual (two gametes) later. Every read is
classified into exactly one imputation outcome (these sum to 100%):

  recall          : placed at the correct locus AND the true founder is in the set.
  founder-dropout : placed at the correct locus but the true founder is MISSING.
  misplaced       : placed at a locus overlapping none of the read's true segments
                    (evidence emitted at the wrong place).
  unmapped        : no confident placement (UNPLACED/MULTI) -- evidence lost.

**founder-dropout and misplaced are the errors that hurt imputation** (the true
founder loses support, or a wrong locus gains false support). The rest is a
*secondary* specificity measure, not a founder-path error.

For recall reads, the non-source founders in the emitted set are split (the true
founder is already present, so its recovery is unaffected):
  IBS founders : non-source founders whose sequence is IDENTICAL over the read
                 span (identity-by-state) -- inherent ambiguity, not an error.
  spurious     : founders added only because sequencing error shortened the exact
                 match so it also matched them -- not truly consistent with the
                 full read. A specificity cost the CRF's multihot largely absorbs;
                 it does NOT drop the true founder.
  missed-IBS   : IBS founders we failed to emit (with strict attribution, ~0).
  resolution / exact : Jaccard(emitted, IBS+source) and fraction where equal.

Placement uses interval overlap (not exact base): refmap extrapolates a full-read
[cL,cR) from the exact core, so an error read's coordinate can sit inside (or, on
the '-' strand, just beside) its true span. Overlap is the robust criterion.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "sim"))
import simlib


def parse_gametes_header(ps4g_path):
    """Return (name2idx, idx2name, pos_sets) from a PS4G file.
    pos_sets[(contig,pos)] = list of frozenset(gamete indices)."""
    name2idx, idx2name = {}, {}
    pos_sets = {}
    in_data = False
    with open(ps4g_path) as f:
        for line in f:
            if line.startswith("#"):
                # gamete rows look like  #B73\t0\t4969  (name may carry :sub)
                body = line[1:].rstrip("\n")
                parts = body.split("\t")
                if len(parts) == 3 and parts[1].isdigit():
                    name = parts[0].split(":")[0]
                    idx = int(parts[1])
                    name2idx[name] = idx
                    idx2name[idx] = name
                continue
            if line.startswith("gameteSet"):
                in_data = True
                continue
            if not in_data:
                continue
            gs, contig, pos, cnt = line.rstrip("\n").split("\t")
            s = frozenset(int(x) for x in gs.split(","))
            pos_sets.setdefault((contig, int(pos)), []).append(s)
    return name2idx, idx2name, pos_sets


def parse_refmap_table(path, name2idx):
    """qname -> dict(status, contig, coord, carriers:set(idx))."""
    out = {}
    with open(path) as f:
        for line in f:
            F = line.rstrip("\n").split("\t")
            if len(F) < 10 or F[0] == "qname":
                continue
            qname, qlen, status, ncar, carriers, refname, strand, cL, cR, span = F[:10]
            contig = simlib.contig_of(refname) if refname != "." else None
            coord = int(cL) if cL not in (".", "") else None
            end = int(cR) if cR not in (".", "") else None
            cset = set()
            if carriers != ".":
                for tok in carriers.split(","):
                    nm = tok.split(":")[0]
                    if nm in name2idx:
                        cset.add(name2idx[nm])
            out[qname] = {"status": status, "contig": contig, "coord": coord,
                          "end": end, "carriers": cset}
    return out


def parse_truth(path, name2idx):
    """qname -> dict(source_idx, segments:[(contig,lo,hi)], oracle:set, has_error)."""
    out = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            F = line.rstrip("\n").split("\t")
            (rid, src, tx, strand, segstr, njunc, errpos, oraclestr, is_pav) = F
            segments = []
            for seg in segstr.split(";"):
                contig, span = seg.rsplit(":", 1)
                lo, hi = span.split("-")
                segments.append((contig, int(lo), int(hi)))
            # whole-read oracle = INTERSECTION of the per-segment oracle sets: a
            # spliced read came from one founder that must match every exon it
            # spans, so the founders consistent with the whole read are those in
            # all segments (within-exon reads have one segment -> unchanged). This
            # is the reference the chaining emitter's intersected set is scored
            # against; it also correctly flags a stock exon-only set that includes
            # a founder present in one exon but not the others as spurious.
            seg_oracles = [set(int(x) for x in part.split(",") if x != "")
                           for part in oraclestr.split("|")]
            oracle = set.intersection(*seg_oracles) if seg_oracles else set()
            out[rid] = {"source_idx": name2idx.get(src),
                        "segments": segments,
                        "oracle": oracle,
                        "seg_oracles": seg_oracles,
                        "has_error": errpos != "",
                        "is_pav": is_pav == "1"}
    return out


def parse_ps4g_per_read(path):
    """The per-read PS4G file (refmap --ps4g-per-read, or the chaining emitter):
    readName -> list of (contig, pos, frozenset(gametes)). A read may emit more
    than one row (one per exon segment, from the chaining emitter); stock refmap
    emits one. This is the exact set a read contributes, so attribution is strict."""
    out = {}
    with open(path) as f:
        f.readline()  # header
        for line in f:
            rid, contig, pos, gs = line.rstrip("\n").split("\t")
            s = frozenset(int(x) for x in gs.split(",")) if gs else frozenset()
            out.setdefault(rid, []).append((contig, int(pos), s))
    return out


def placed_at_locus(rec, segments):
    """True when the read's emitted reference interval overlaps one of its true
    exon segments on the matching contig. refmap extrapolates a full-read
    [cL,cR) interval from the exact core, so for error-containing reads (esp. on
    the '-' strand) the interval is shifted by up to the read length but still
    overlaps the true locus; overlap is the robust, tight criterion."""
    contig, lo0 = rec["contig"], rec["coord"]
    if contig is None or lo0 is None:
        return False
    hi0 = rec["end"] if rec["end"] is not None else lo0 + 1
    for (c, lo, hi) in segments:
        if c == contig and lo0 < hi and lo < hi0:
            return True
    return False


def jaccard(a, b):
    if not a and not b:
        return 1.0
    u = a | b
    return len(a & b) / len(u) if u else 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--ps4g-per-read", dest="ps4g_per_read", required=True,
                    help="per-read PS4G file (refmap --ps4g-per-read or the chaining "
                         "emitter): exact per-read gameteSet(s), for strict attribution")
    ap.add_argument("--gametes", help="gametes.tsv (gamete index map); used for names "
                    "when --ps4g is absent, else validated against the PS4G header")
    ap.add_argument("--ps4g", help="aggregated PS4G (stock refmap) for header + row "
                    "validation; optional for the chaining emitter")
    ap.add_argument("--refmap", help="stock refmap per-read table for exact [cL,cR] "
                    "placement geometry; optional (else reconstructed from pos+len)")
    a = ap.parse_args()

    pos_sets = {}
    if a.ps4g:
        name2idx, idx2name, pos_sets = parse_gametes_header(a.ps4g)
    elif a.gametes:
        name2idx = {}
        for line in open(a.gametes):
            idx, name = line.rstrip("\n").split("\t")
            name2idx[name] = int(idx)
    else:
        ap.error("need --ps4g or --gametes for the gamete index map")
    table = parse_refmap_table(a.refmap, name2idx) if a.refmap else {}
    truth = parse_truth(a.truth, name2idx)
    per_read = parse_ps4g_per_read(a.ps4g_per_read)  # readName -> [(contig, pos, set), ...]

    # the individual's founder set: reads are sampled from one gamete now
    # (homozygous), a heterozygous individual (two gametes) later. Non-founder
    # emissions are judged against this set; = {B73} for a single-source run.
    founders = {t["source_idx"] for t in truth.values() if t["source_idx"] is not None}

    # --- PS4G validation
    problems = []
    if a.gametes:
        want = {}
        for line in open(a.gametes):
            idx, name = line.rstrip("\n").split("\t")
            want[name] = int(idx)
        if want != name2idx:
            problems.append("gamete index map mismatch: %s vs PS4G %s"
                            % (want, name2idx))
    # every PS4G row position should fall within (or within a read length of, to
    # allow refmap's core->full-read extrapolation shift) some true read span;
    # rows far from any true locus would be spurious (e.g. intronic) evidence.
    truth_spans = {}
    pad = 0
    for t in truth.values():
        for (c, lo, hi) in t["segments"]:
            truth_spans.setdefault(c, []).append((lo, hi))
            pad = max(pad, hi - lo)
    rows_total = rows_in_span = 0
    for (contig, pos), sets in pos_sets.items():
        for _ in sets:
            rows_total += 1
            if any(lo - pad <= pos < hi + pad for (lo, hi) in truth_spans.get(contig, [])):
                rows_in_span += 1

    # --- per-read scoring. The emitted gamete set is the *exact* set the read
    # contributed, read from the per-read PS4G file (strict attribution) -- no
    # lookup into the aggregated PS4G, so co-located reads never mix.
    def bucket():
        return {"n": 0,
                # placement partition (mutually exclusive, sum == n)
                "recall": 0, "dropout": 0, "wrong_region": 0, "unmapped": 0,
                # set quality, accumulated over recall reads only
                "res_sum": 0.0, "exact": 0,
                "over_reads": 0, "over_sum": 0,     # algorithmic false founders (vs oracle)
                "under_reads": 0, "under_sum": 0,   # missed truly-consistent founders (vs oracle)
                "ambig_sum": 0}                     # inherent identical-gamete ambiguity (vs source)

    allb, clean, noisy = bucket(), bucket(), bucket()
    for rid, t in truth.items():
        rows = per_read.get(rid)    # [(contig, pos, set), ...] or None
        rec = table.get(rid)
        b = clean if not t["has_error"] else noisy
        for bb in (allb, b):
            bb["n"] += 1
        # --- classify placement into exactly one bucket
        if not rows:                # not EXACT/PLACED -> no evidence emitted
            for bb in (allb, b):
                bb["unmapped"] += 1
            continue
        es = set().union(*[set(r[2]) for r in rows])   # emitted set (chain rows share it)
        # placement: prefer the stock table's exact [cL,cR); else reconstruct
        # [pos, pos+read_len) and accept if ANY emitted row overlaps a true segment.
        if rec is not None and rec["coord"] is not None:
            at_locus = placed_at_locus(rec, t["segments"])
        else:
            L = sum(hi - lo for _, lo, hi in t["segments"]) or 1
            at_locus = any(rc == c and rp < hi and lo < rp + L
                           for (rc, rp, _) in rows for (c, lo, hi) in t["segments"])
        if not at_locus:
            for bb in (allb, b):
                bb["wrong_region"] += 1
            continue
        if t["source_idx"] not in es:
            for bb in (allb, b):
                bb["dropout"] += 1
            continue
        # --- recall read: measure set quality (source excluded, both references)
        O = t["oracle"]
        e_nonf = es - founders                 # non-founder emissions
        over = e_nonf - O                      # false founders (not truly consistent)
        ambig = e_nonf & O                     # truly consistent -> inherent ambiguity
        under = (O - founders) - es            # truly-consistent founders we missed
        res = jaccard(es, O)
        for bb in (allb, b):
            bb["recall"] += 1
            bb["res_sum"] += res
            if es == O:
                bb["exact"] += 1
            if over:
                bb["over_reads"] += 1
            bb["over_sum"] += len(over)
            if under:
                bb["under_reads"] += 1
            bb["under_sum"] += len(under)
            bb["ambig_sum"] += len(ambig)

    def report(name, b):
        if b["n"] == 0:
            return
        n = b["n"]
        pct = lambda k: 100.0 * b[k] / n
        print("  %-11s n=%-6d recall=%5.1f%%  founder-dropout=%4.1f%%  misplaced=%4.1f%%  unmapped=%4.1f%%"
              % (name, n, pct("recall"), pct("dropout"), pct("wrong_region"), pct("unmapped")))
        r = b["recall"]
        if r:
            per = lambda k: b[k] / r
            rate = lambda k: 100.0 * b[k] / r
            print("  %-11s   IBS=%.3f/rd  spurious=%4.1f%%(%.3f/rd)  missed-IBS=%4.1f%%(%.3f/rd)"
                  "  resolution=%.3f exact=%4.1f%%"
                  % ("", per("ambig_sum"),
                     rate("over_reads"), per("over_sum"),
                     rate("under_reads"), per("under_sum"),
                     b["res_sum"] / r, 100.0 * b["exact"] / r))

    print("== imputation outcome (per read; recall/founder-dropout/misplaced/unmapped sum to 100%) ==")
    print("reads=%d  emitted (per-read PS4G file)=%d" % (len(truth), len(per_read)))
    report("ALL", allb)
    report("error-free", clean)
    report("with-error", noisy)
    print("  CRITICAL for imputation: founder-dropout (true founder lost at the right locus)")
    print("  and misplaced (evidence at the wrong locus). The 2nd line is SECONDARY founder-")
    print("  set breadth among recovered reads (true founder already present): IBS = non-")
    print("  source founders identical over the span (inherent); spurious = founders added by")
    print("  error, not truly consistent -- a specificity cost the CRF's multihot absorbs.")
    if a.ps4g:
        print("== PS4G validation ==")
        print("  gamete index map: %s"
              % ("OK" if not problems else "; ".join(problems)))
        print("  rows within a true read span: %d/%d (%.1f%%)"
              % (rows_in_span, rows_total,
                 100.0 * rows_in_span / rows_total if rows_total else 0.0))
    if problems:
        sys.exit(1)


if __name__ == "__main__":
    main()
