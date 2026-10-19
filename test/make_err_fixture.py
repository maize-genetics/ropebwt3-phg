#!/usr/bin/env python3
"""Generate the --err-max-events integration reads in docs/examples/err/ (committed; rerun only to change it).

Uses the --paired fixture pangenome (docs/examples/pe/pangenome.fa; B73_chr1 = 20 kb random R,
Oh43_chr1 = R unchanged, Mo17/Ki3 = R plus insertions). Single-end reads are 150 bp copies of
R[3000:3150] etc. (a region every genome carries unchanged), edited by hand so the expected
event count is known:
  E0  exact                                   -> 0 events
  E1  one substitution at read offset 75      -> 1
  E2  substitutions at 40 and 110             -> 2
  E3  substitutions at 30, 75 and 120         -> 3
  I1  one extra base inserted at offset 80    -> 1 (read insertion, not 2+ substitutions)
  D1  one base deleted at offset 80           -> 1 (read deletion)
  L1  substitution at offset 3 (left end)     -> 1
  R1  substitution at offset 147 (right end)  -> 1
Pairs (r1/r2, mate 2 reverse-complemented, 400 bp fragments in R[4000:9000]):
  Q0 0+0, Q1 1+0, Q2 1+1, Q3 2+1 events.
"""
import os

here = os.path.dirname(os.path.abspath(__file__))
fa = os.path.join(here, "..", "docs", "examples", "pe", "pangenome.fa")
out = os.path.join(here, "..", "docs", "examples", "err")
seqs, name = {}, None
for line in open(fa):
    line = line.strip()
    if line.startswith(">"): name = line[1:]; seqs[name] = []
    else: seqs[name].append(line)
R = "".join(seqs["B73_chr1"])
comp = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(comp)[::-1]
def sub(s, offs):
    s = list(s)
    for o in offs: s[o] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[o]]
    return "".join(s)

base = lambda st: R[st:st + 150]
reads = [("E0", base(3000)), ("E1", sub(base(3200), [75])), ("E2", sub(base(3400), [40, 110])),
         ("E3", sub(base(3600), [30, 75, 120])),
         ("I1", base(3800)[:80] + {"A": "C", "C": "G", "G": "T", "T": "A"}[R[3880]] + base(3800)[80:149]),
         ("D1", R[4000:4080] + R[4081:4151]), ("L1", sub(base(4200), [3])), ("R1", sub(base(4400), [147]))]
with open(os.path.join(out, "se.fa"), "w") as f:
    for n, s in reads: f.write(f">{n}\n{s}\n")
pairs = [("Q0", [], []), ("Q1", [60], []), ("Q2", [60], [50]), ("Q3", [30, 110], [50])]
with open(os.path.join(out, "r1.fa"), "w") as f1, open(os.path.join(out, "r2.fa"), "w") as f2:
    for i, (n, e1, e2) in enumerate(pairs):
        st = 5000 + 900 * i
        m1, m2 = R[st:st + 150], rc(R[st + 250:st + 400])
        f1.write(f">{n}/1\n{sub(m1, e1)}\n"); f2.write(f">{n}/2\n{sub(m2, e2)}\n")
