#!/usr/bin/env python3
"""Generate the --paired integration fixture in docs/examples/pe/ (committed; rerun only to change it).

Pangenome (4 genomes, one contig each):
  B73_chr1   R: 20,000 random bp (the reference); R[11000:11300] is a copy of R[19000:19300]
  Mo17_chr1  R with a 1,000 bp insertion I at R position 10,000 and a 400 bp element T at R 9,500
  Ki3_chr1   R with I at 10,000, SNPs at R positions 10,320 / 10,380 / 10,440 (right of the
             insertion), and the element T + a Ki3-only element U at R 15,000 (a TE copy elsewhere)
  Oh43_chr1  R, unchanged
Gamete indices (B73 first, then alphabetical): B73=0 Ki3=1 Mo17=2 Oh43=3.

Pairs (mate 1 forward, mate 2 reverse-complemented, 100 bp unless noted):
  iso00..iso39  fragments of 380..419 bp from R[1000:9000] -> 40 EXACT/EXACT pairs, spans 380..419
                (insert-size estimate: mean 399.5)
  P1   mate 1 = I[400:550] (Mo17+Ki3 only -> PLACED via lift), mate 2 = R[10300:10450] (B73, Mo17,
       Oh43 -- Ki3 carries SNPs there -> EXACT). One CONC event, founder set = {Mo17} = "2".
  P2   both mates inside R[12000:12500] (all four genomes) -> CONC, set {0,1,2,3}
  P3   mate 1 = R[14000:14100], mate 2 = random junk -> SOLO, set {0,1,2,3}
  P4   mate 1 = R[1000:1100], mate 2 = R[18000:18100] -> DISC (17 kb apart): two events
  P5   mate 1 = T[100:250] (Mo17 projects ~9,200, Ki3 ~16,100; the upper median picks Ki3's),
       mate 2 = R[9600:9700] (EXACT) -> mate 1 MOVED next to mate 2 (RESC); Ki3's copy is
       elsewhere, so the event's founder set is {Mo17} = "2"
  P6   mate 1 = R[19000:19100] (two reference copies: 11000 and 19000), mate 2 = R[11350:11450]
       -> mate 1 ends at the copy next to its mate (11000): CONC, or RESC/MOVED if refmap reported
       the other copy first
  P7   mate 1 = U[50:200] (Ki3 only, projects ~16,450), mate 2 = R[6000:6100] (EXACT) -> DISC;
       with --pe-trust-exact mate 1 is MOVED onto mate 2 (6000)
"""
import os, random

rng = random.Random(20261014)
comp = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(comp)[::-1]
def rand_seq(n): return "".join(rng.choice("ACGT") for _ in range(n))

R = rand_seq(20000)
I = rand_seq(1000)
T = rand_seq(400)
U = rand_seq(300)
R = R[:11000] + R[19000:19300] + R[11300:]  # segmental duplicate inside the reference, 8 kb from its source
mo17 = R[:9500] + T + R[9500:10000] + I + R[10000:]
ki3 = list(R[:10000] + I + R[10000:])
for rp in (10320, 10380, 10440):
    q = rp + 1000  # position in the assembly (after the insertion)
    ki3[q] = {"A": "C", "C": "G", "G": "T", "T": "A"}[ki3[q]]
ki3 = "".join(ki3)
ki3 = ki3[:15000 + 1000] + T + U + ki3[15000 + 1000:]

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "examples", "pe")
os.makedirs(out, exist_ok=True)

def wrap(s, w=80): return "\n".join(s[i:i + w] for i in range(0, len(s), w))
with open(os.path.join(out, "pangenome.fa"), "w") as f:
    for name, s in (("B73_chr1", R), ("Mo17_chr1", mo17), ("Ki3_chr1", ki3), ("Oh43_chr1", R)):
        f.write(f">{name}\n{wrap(s)}\n")

pairs = []
for i in range(40):
    L = 380 + (i * 13) % 40          # all of 380..419, shuffled
    s = 1000 + i * 190
    pairs.append((f"iso{i:02d}", R[s:s + 100], rc(R[s + L - 100:s + L])))
pairs.append(("P1", I[400:550], rc(R[10300:10450])))
pairs.append(("P2", R[12000:12100], rc(R[12300:12400])))
pairs.append(("P3", R[14000:14100], rand_seq(100)))
pairs.append(("P4", R[1000:1100], rc(R[18000:18100])))
pairs.append(("P5", T[100:250], rc(R[9600:9700])))
pairs.append(("P6", R[19000:19100], rc(R[11350:11450])))
pairs.append(("P7", U[50:200], rc(R[6000:6100])))

with open(os.path.join(out, "r1.fa"), "w") as f1, open(os.path.join(out, "r2.fa"), "w") as f2, \
     open(os.path.join(out, "interleaved.fa"), "w") as fi:
    for n, a, b in pairs:
        f1.write(f">{n}/1\n{a}\n"); f2.write(f">{n}/2\n{b}\n")
        fi.write(f">{n}/1\n{a}\n>{n}/2\n{b}\n")
print("wrote", out)
