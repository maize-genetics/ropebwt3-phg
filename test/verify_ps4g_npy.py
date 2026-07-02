#!/usr/bin/env python3
"""Verify that a refmap --npy array corresponds exactly to its --ps4g file.

Usage: verify_ps4g_npy.py <file.ps4g> <file.npy>
  (expects <file.npy>.bins.tsv and <file.npy>.gametes.tsv alongside <file.npy>)

Checks, in order:
  1. row counts match: PS4G data rows == npy rows == bins.tsv rows
  2. gamete counts match: PS4G #gamete rows == gametes.tsv rows == npy cols - 2
  3. gamete name/index tables are identical between the PS4G header and gametes.tsv
  4. per row: bins.tsv (contig, bin) == the PS4G row's (refContig, refPosBinned)
  5. per row: the set of nonzero gamete columns in npy == the PS4G row's gameteSet
     (this is the property that matters: a row must never be aggregated across
     gameteSets, and every gamete the npy row claims support for must be one
     PS4G actually listed, and vice versa)
  6. per row: nonzero cell values equal the PS4G row's count (count mode) or 1
     (binary mode, auto-detected)
  7. per gamete: the PS4G header's total count equals the sum, over rows
     containing that gamete, of the row's count (an internal PS4G invariant)

Exits 0 and prints "OK" with a short summary if everything matches, else exits
1 and prints exactly what mismatched.
"""
import sys
import numpy as np


def parse_ps4g(path):
    gamete_name = {}
    gamete_hdr_count = {}
    rows = []  # (gamete_set: sorted tuple[int], contig: str, bin: int, count: int)
    with open(path) as f:
        in_data = False
        for line in f:
            line = line.rstrip("\n")
            if not in_data:
                if line == "gameteSet\trefContig\trefPosBinned\tcount":
                    in_data = True
                    continue
                if line.startswith("#") and line not in ("#PS4G",) and not line.startswith("#version=") \
                        and not line.startswith("#Command:") and not line.startswith("#TotalUniqueCounts:") \
                        and line != "#gamete\tgameteIndex\tcount":
                    name, idx, cnt = line[1:].split("\t")
                    gamete_name[int(idx)] = name
                    gamete_hdr_count[int(idx)] = int(cnt)
            else:
                gset, contig, binpos, count = line.split("\t")
                gamete_set = tuple(sorted(int(x) for x in gset.split(",")))
                rows.append((gamete_set, contig, int(binpos), int(count)))
    return gamete_name, gamete_hdr_count, rows


def main():
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} <file.ps4g> <file.npy>", file=sys.stderr)
        return 2
    ps4g_fn, npy_fn = sys.argv[1], sys.argv[2]

    gamete_name, gamete_hdr_count, rows = parse_ps4g(ps4g_fn)
    n_gamete_ps4g = len(gamete_name)
    R = len(rows)

    mat = np.load(npy_fn)
    bins_tsv = [l.rstrip("\n").split("\t") for l in open(npy_fn + ".bins.tsv")][1:]
    gametes_tsv = [l.rstrip("\n").split("\t") for l in open(npy_fn + ".gametes.tsv")][1:]

    errors = []

    # 1. row counts
    if not (R == mat.shape[0] == len(bins_tsv)):
        errors.append(f"row count mismatch: ps4g={R} npy={mat.shape[0]} bins.tsv={len(bins_tsv)}")

    # 2/3. gamete tables
    G = mat.shape[1] - 2
    if not (n_gamete_ps4g == len(gametes_tsv) == G):
        errors.append(f"gamete count mismatch: ps4g header={n_gamete_ps4g} gametes.tsv={len(gametes_tsv)} npy_cols-2={G}")
    gametes_tsv_map = {int(idx): name for idx, name in gametes_tsv}
    if gametes_tsv_map != gamete_name:
        errors.append(f"gamete name/index table differs between ps4g header and gametes.tsv: {gamete_name} vs {gametes_tsv_map}")

    if errors:
        print("MISMATCH:")
        for e in errors:
            print(f"  - {e}")
        return 1

    # 4. per-row bins.tsv vs ps4g contig/bin
    bad_bins = []
    for i, (gset, contig, binpos, count) in enumerate(rows):
        bt_contig, bt_bin = bins_tsv[i][1], int(bins_tsv[i][2])
        if bt_contig != contig or bt_bin != binpos:
            bad_bins.append((i, (contig, binpos), (bt_contig, bt_bin)))
    if bad_bins:
        print(f"MISMATCH: {len(bad_bins)} row(s) where bins.tsv (contig,bin) != ps4g (refContig,refPosBinned)")
        for i, want, got in bad_bins[:10]:
            print(f"  - row {i}: ps4g={want} bins.tsv={got}")
        return 1

    # 5/6. nonzero gamete columns == gameteSet, and values match count or binary(1)
    row_idx = np.repeat(np.arange(R), [len(gset) for gset, *_ in rows])
    col_idx = np.concatenate([np.array(gset, dtype=np.int64) for gset, *_ in rows]) if R else np.array([], dtype=np.int64)
    ps4g_mask = np.zeros((R, G), dtype=bool)
    ps4g_mask[row_idx, col_idx] = True

    npy_gamete_cols = mat[:, :G]
    npy_mask = npy_gamete_cols != 0

    if not np.array_equal(ps4g_mask, npy_mask):
        mismatched_rows = np.nonzero(np.any(ps4g_mask != npy_mask, axis=1))[0]
        print(f"MISMATCH: {len(mismatched_rows)} row(s) where the npy nonzero-gamete set != the ps4g gameteSet")
        for i in mismatched_rows[:10]:
            npy_set = tuple(int(x) for x in np.nonzero(npy_mask[i])[0])
            print(f"  - row {i}: ps4g gameteSet={rows[i][0]} npy nonzero cols={npy_set}")
        return 1

    counts = np.array([count for *_, count in rows], dtype=np.int64)
    nonzero_vals = npy_gamete_cols[ps4g_mask]
    is_binary = bool(R) and np.all(nonzero_vals == 1) and not np.all(counts == 1)
    if is_binary:
        expected = np.ones_like(npy_gamete_cols)
    else:
        expected = np.broadcast_to(counts[:, None], npy_gamete_cols.shape)
    actual_masked = np.where(ps4g_mask, npy_gamete_cols, 0)
    expected_masked = np.where(ps4g_mask, expected, 0)
    if not np.array_equal(actual_masked, expected_masked):
        bad = np.nonzero(np.any(actual_masked != expected_masked, axis=1))[0]
        print(f"MISMATCH: {len(bad)} row(s) where npy cell values don't match ps4g count ({'binary' if is_binary else 'count'} mode)")
        for i in bad[:10]:
            print(f"  - row {i}: ps4g count={rows[i][3]} npy row={npy_gamete_cols[i].tolist()}")
        return 1

    # 7. per-gamete header total == sum of row counts containing that gamete
    gamete_sum = np.zeros(G, dtype=np.int64)
    np.add.at(gamete_sum, col_idx, np.repeat(counts, [len(gset) for gset, *_ in rows]))
    bad_totals = []
    for idx, hdr_count in gamete_hdr_count.items():
        if gamete_sum[idx] != hdr_count:
            bad_totals.append((idx, gamete_name[idx], hdr_count, int(gamete_sum[idx])))
    if bad_totals:
        print(f"MISMATCH: {len(bad_totals)} gamete(s) where the #gamete header count != sum of row counts")
        for idx, name, hdr, actual in bad_totals[:10]:
            print(f"  - gamete {idx} ({name}): header={hdr} actual_sum={actual}")
        return 1

    print(f"OK: {R} rows, {G} gametes, {'binary' if is_binary else 'count'} mode -- "
          f"npy corresponds exactly to {ps4g_fn}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
