#!/usr/bin/env python3
"""Paired-end workflow diagram (docs/paired-end/plan_diagram.png). Regenerate whenever the
stage list in PLAN.md changes:  python3 docs/paired-end/plan_diagram.py  (needs rsvg-convert)."""
import os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
COL = {"reused": ("#eceff1", "#90a4ae"), "new": ("#e3f6e5", "#43a047"), "mod": ("#efe6fb", "#7e57c2"),
       "later": ("#ffffff", "#b0bec5")}
TAG = {"reused": "REUSED", "new": "NEW", "mod": "MODIFIED", "later": "LATER PHASE"}

STAGES = [  # (title, detail, kind)
    ("Input: --paired R1 R2 | interleaved", "main_search", "mod"),
    ("Batch reader: whole pairs, mate-name check", "worker_pipeline step 0", "mod"),
    ("Per-read placement (unchanged per read)", "refmap_query", "reused"),
    ("Keep alternatives: lift projections, EXACT copies", "refmap_place_lift, refmap_query", "mod"),
    ("Insert size: per batch, or --pe-isize", "rb3_pe_isize_estimate", "new"),
    ("Classify pair: CONC RESC SOLO DISC NONE", "rb3_pe_resolve", "new"),
    ("Rescue / tie-break: move a mate next to the other", "rb3_pe_resolve, rb3_pe_near_gametes", "new"),
    ("stdout: + pair-class column", "write_refmap", "mod"),
    ("1 event per pair, founders = intersection", "refmap_pair_event", "mod"),
    ("PS4G / npy finalize (layout unchanged)", "rb3_ps4g_npy_finalize", "reused"),
]
LATER = [("Phase 2: chain", "link SMEM chains across mates"),
         ("Phase 3: mem / sw", "bwa-style mate rescue + pair score"),
         ("Simulator", "paired training reads (separate)")]

W, BX, BW, BH, GAP, TOP = 980, 40, 560, 54, 18, 90
H = TOP + len(STAGES) * (BH + GAP) + 150


def box(x, y, w, h, kind, title, detail, dashed=False):
    fill, stroke = COL[kind]
    dash = ' stroke-dasharray="7,5"' if dashed else ""
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="9" fill="{fill}" stroke="{stroke}" stroke-width="2"{dash}/>'
            f'<text x="{x+14}" y="{y+23}" font-family="Liberation Sans" font-size="16" font-weight="bold" fill="#263238">{title}</text>'
            f'<text x="{x+14}" y="{y+43}" font-family="Liberation Mono" font-size="13" fill="#546e7a">{detail}</text>')


svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
       '<defs><marker id="a" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto">'
       '<path d="M0,0 L10,4 L0,8 z" fill="#78909c"/></marker></defs>',
       '<text x="40" y="44" font-family="Liberation Sans" font-size="24" font-weight="bold" fill="#263238">'
       'Paired-end reads in ropebwt3 refmap (phase 1)</text>',
       '<text x="40" y="70" font-family="Liberation Sans" font-size="14" fill="#546e7a">'
       'Opt-in with --paired; without it every stage below runs exactly as before.</text>']
for i, (t, d, k) in enumerate(STAGES):
    y = TOP + i * (BH + GAP)
    svg.append(box(BX, y, BW, BH, k, t, d))
    fill, stroke = COL[k]
    svg.append(f'<text x="{BX+BW-12}" y="{y+23}" text-anchor="end" font-family="Liberation Sans" font-size="11" '
               f'font-weight="bold" fill="{stroke}">{TAG[k]}</text>')
    if i + 1 < len(STAGES):
        svg.append(f'<line x1="{BX+BW/2}" y1="{y+BH}" x2="{BX+BW/2}" y2="{y+BH+GAP-2}" stroke="#78909c" stroke-width="2" marker-end="url(#a)"/>')

# later phases
LX, LW = 640, 300
svg.append(f'<text x="{LX}" y="{TOP-10}" font-family="Liberation Sans" font-size="14" font-weight="bold" fill="#78909c">Not in phase 1</text>')
for i, (t, d) in enumerate(LATER):
    y = TOP + i * (BH + GAP)
    svg.append(box(LX, y, LW, BH, "later", t, d, dashed=True))

# legend
ly = TOP + 4 * (BH + GAP) + 10
for i, k in enumerate(["reused", "new", "mod"]):
    fill, stroke = COL[k]
    y = ly + i * 30
    svg.append(f'<rect x="{LX}" y="{y}" width="22" height="18" rx="4" fill="{fill}" stroke="{stroke}" stroke-width="2"/>'
               f'<text x="{LX+32}" y="{y+14}" font-family="Liberation Sans" font-size="14" fill="#263238">{TAG[k].capitalize()}</text>')

# acceptance gate
gy = TOP + len(STAGES) * (BH + GAP) + 6
svg.append(f'<rect x="{BX}" y="{gy}" width="{W-2*BX}" height="96" rx="9" fill="#fff8e1" stroke="#f9a825" stroke-width="2"/>')
for i, line in enumerate(["Acceptance",
                          "Single-end output byte-identical to the parent binaries (refmap + chain, real reads)",
                          "Unit tests (test_pairing) + fixture integration (run_integration section 11, valgrind)",
                          "Maize wgsim pairs (Oh43 in-panel, Ia453 held-out): placement + event scoring, SE vs PE"]):
    bold = 'font-weight="bold" ' if i == 0 else ""
    svg.append(f'<text x="{BX+16}" y="{gy+24+i*21}" font-family="Liberation Sans" font-size="{15 if i == 0 else 13}" '
               f'{bold}fill="#5d4037">{line}</text>')
svg.append("</svg>")

svg_fn = os.path.join(HERE, "plan_diagram.svg")
open(svg_fn, "w").write("\n".join(svg))
subprocess.run(["rsvg-convert", "-o", os.path.join(HERE, "plan_diagram.png"), svg_fn], check=True)
os.remove(svg_fn)
print("wrote", os.path.join(HERE, "plan_diagram.png"))
