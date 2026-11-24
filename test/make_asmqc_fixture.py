#!/usr/bin/env python3
"""Synthetic reference/assembly pair with planted events at known positions, for test/run_asmqc.sh.

Writes <out>/REF.fa (taxon REF), <out>/CAND.fa (taxon CAND) and <out>/truth.tsv. Deterministic (fixed seed),
standard library only, ~2 Mb of sequence.

REF chr1 (1.2 Mb) vs CAND chr1:
  SNPs at 500 + 997*i in [0, 400 kb)            -> exact SNP positions
  19 small indels (3 bp ins / 3 bp del, alternating) at 410, 420, ..., 590 kb
  a 3 kb replacement at 700,000                  -> one replacement (3000 ref / 3000 asm bp)
  a 5 kb insertion at 800,000
  a 20 kb deletion at 900,000
  an inversion of [1,000,000, 1,100,000)
REF chr2 (600 kb): its first 100 kb is moved to the end of CAND chr1 (translocation)
REF chr3 (200 kb): carries a 30 kb LTR element at 60,000 (LTR-internal-LTR); CAND chr3 has only a solo LTR
  (10% diverged) there                          -> ~150 SNPs + one ~28.5 kb deletion
"""
import os, random, sys

out = sys.argv[1] if len(sys.argv) > 1 else '.'
os.makedirs(out, exist_ok=True)
rng = random.Random(20261124)
def rnd(n): return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(s): return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
def sub(c): return rng.choice([b for b in 'ACGT' if b != c])

c1, c2, c3 = rnd(1_200_000), rnd(600_000), rnd(170_000)
ltr, internal = rnd(1500), rnd(27000)
te = ltr + internal + ''.join(sub(b) if rng.random() < 0.02 else b for b in ltr)
c3 = c3[:60000] + te + c3[60000:]
truth = []

a = list(c1)
for p in range(500, 400_000, 997):
    a[p] = sub(a[p]); truth.append(('SNP', 'chr1', p + 1))
a = ''.join(a)
pieces, last = [], 0
for i, p in enumerate(range(410_000, 600_000, 10_000)):   
    pieces.append(a[last:p])
    if i % 2 == 0: pieces.append(rnd(3)); truth.append(('INS3', 'chr1', p + 1)); last = p
    else: truth.append(('DEL3', 'chr1', p + 1)); last = p + 3
pieces.append(a[last:])
a = ''.join(pieces)
shift = sum(3 if t[0] == 'INS3' else -3 for t in truth if t[0] in ('INS3', 'DEL3'))  # net length change before 700 kb
s700 = 700_000 + shift
a = a[:s700] + rnd(3000) + a[s700 + 3000:]; truth.append(('REP3000', 'chr1', 700_001))
s800 = 800_000 + shift
a = a[:s800] + rnd(5000) + a[s800:]; truth.append(('INS5000', 'chr1', 800_001))
s900 = 900_000 + shift + 5000
a = a[:s900] + a[s900 + 20000:]; truth.append(('DEL20000', 'chr1', 900_001))
s1m = 1_000_000 + shift + 5000 - 20000
a = a[:s1m] + rc(a[s1m:s1m + 100_000]) + a[s1m + 100_000:]; truth.append(('INV100000', 'chr1', 1_000_001))
cand1 = a + c2[:100_000]; truth.append(('TRANSLOC100000', 'chr2', 1))
cand2 = c2[100_000:]
solo = ''.join(sub(b) if rng.random() < 0.10 else b for b in ltr)
cand3 = c3[:60000] + solo + c3[60000 + len(te):]; truth.append(('SOLO_LTR_DEL', 'chr3', 60_001))

def wfa(fn, recs):
    with open(fn, 'w') as f:
        for name, s in recs:
            f.write(f'>{name}\n')
            for i in range(0, len(s), 100): f.write(s[i:i + 100] + '\n')
wfa(os.path.join(out, 'REF.fa'), [('chr1', c1), ('chr2', c2), ('chr3', c3)])
wfa(os.path.join(out, 'CAND.fa'), [('CAND_chr1', cand1), ('CAND_chr2', cand2), ('CAND_chr3', cand3)])
with open(os.path.join(out, 'truth.tsv'), 'w') as f:
    for t in truth: f.write('\t'.join(map(str, t)) + '\n')
