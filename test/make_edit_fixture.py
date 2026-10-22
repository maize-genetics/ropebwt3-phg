#!/usr/bin/env python3
"""Generate the edit-mode (--err-noclip / --err-cigar / --min-match-frac) reads in docs/examples/edit/
(committed; rerun only to change it).

Uses the --paired fixture pangenome (docs/examples/pe/pangenome.fa; B73_chr1 = 20 kb random R, Oh43 = R,
Mo17/Ki3 = R plus internal insertions; every genome starts and ends with R's ends). 150 bp reads, hand
edited so the expected edits / '=' bases / extended CIGAR are known exactly. name<TAB>... comments in
docs/examples/edit/expected.tsv: name, legacy events, edits, '=' bases, CIGAR.
"""
import os, random

here = os.path.dirname(os.path.abspath(__file__))
fa = os.path.join(here, "..", "docs", "examples", "pe", "pangenome.fa")
out = os.path.join(here, "..", "docs", "examples", "edit")
os.makedirs(out, exist_ok=True)
seqs, name = {}, None
for line in open(fa):
    line = line.strip()
    if line.startswith(">"): name = line[1:]; seqs[name] = []
    else: seqs[name].append(line)
R = "".join(seqs["B73_chr1"])
N = len(R)
nxt = {"A": "C", "C": "G", "G": "T", "T": "A"}
def sub(s, offs):
    s = list(s)
    for o in offs: s[o] = nxt[s[o]]
    return "".join(s)
def other(c): return nxt[c]
rng = random.Random(7)
def junk(k, avoid):  # k random bases whose first base differs from `avoid` (so the overhang starts with a mismatch)
    s = [rng.choice("ACGT") for _ in range(k)]
    if s[0] == avoid: s[0] = nxt[avoid]
    return "".join(s)

b = lambda st: R[st:st + 150]
reads = []  # (name, seq, legacy_events, edits, eq, cigar)
reads.append(("E0", b(3000), 0, 0, 150, "150="))
reads.append(("S1", sub(b(3200), [75]), 1, 1, 149, "75=1X74="))
reads.append(("S2", sub(b(3400), [40, 110]), 2, 2, 148, "40=1X69=1X39="))
# 1-bp read insertion after read offset 80 (base differs from both neighbours' carrier bases)
st = 3800; ins = [c for c in "ACGT" if c != R[st + 80] and c != R[st + 79]][0]
reads.append(("I1", R[st:st + 80] + ins + R[st + 80:st + 149], 1, 1, 149, "80=1I69="))
# 3-bp read insertion after offset 80
st = 4600; ins3 = "".join(other(R[st + 80 + j]) for j in range(3))
if ins3[0] == R[st + 79]: ins3 = other(ins3[0]) + ins3[1:]
reads.append(("I3", R[st:st + 80] + ins3 + R[st + 80:st + 147], None, 3, 147, "80=3I67="))
# 1-bp and 3-bp read deletions at offset 80
reads.append(("D1", R[4000:4080] + R[4081:4151], 1, 1, 150, "80=1D70="))
st = 5000  # unambiguous placement: the base after the deletion differs from the first deleted one
while R[st + 80] == R[st + 83]: st += 1
reads.append(("D3", R[st:st + 80] + R[st + 83:st + 153], None, 3, 150, "80=3D70="))
# mismatches in the last base at either end
reads.append(("X0", sub(b(5600), [0]), 1, 1, 149, "1X149="))
reads.append(("X149", sub(b(5800), [149]), 1, 1, 149, "149=1X"))
# contig-end overhangs at the right end of every genome (k bp past R's last base), and at the left
for k in (1, 5, 20):
    reads.append((f"OR{k}", R[N - (150 - k):] + junk(k, "A"), 1, k, 150 - k, f"{150 - k}={k}S"))
for k in (1, 5, 20):
    j = junk(k, "A")[::-1]
    reads.append((f"OL{k}", j + R[:150 - k], 1, k, 150 - k, f"{k}S{150 - k}="))
with open(os.path.join(out, "se.fa"), "w") as f:
    for n, s, *_ in reads: f.write(f">{n}\n{s}\n")
with open(os.path.join(out, "expected.tsv"), "w") as f:
    f.write("name\tlegacy_events\tedits\teq\tcigar\n")
    for n, s, le, ed, eq, cg in reads:
        assert len(s) == 150, (n, len(s))
        f.write(f"{n}\t{'.' if le is None else le}\t{ed}\t{eq}\t{cg}\n")
# Pairs (mate 2 reverse-complemented): M0 0+0 edits, M2 1+1 substitutions, M3 2+1, C5 = a clean mate
# plus a mate running 5 bp off the right contig end (pair edits 5, overhang -> --err-clip-reject).
comp = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(comp)[::-1]
pairs = [("M0", R[7000:7150], rc(R[7250:7400])), ("M2", sub(R[7600:7750], [60]), rc(sub(R[7850:8000], [50]))),
         ("M3", sub(R[8200:8350], [30, 110]), rc(sub(R[8450:8600], [50]))),
         ("C5", R[N - 400:N - 250], rc(R[N - 145:] + junk(5, "A")))]
with open(os.path.join(out, "r1.fa"), "w") as f1, open(os.path.join(out, "r2.fa"), "w") as f2:
    for n, m1, m2 in pairs: f1.write(f">{n}/1\n{m1}\n"); f2.write(f">{n}/2\n{m2}\n")
