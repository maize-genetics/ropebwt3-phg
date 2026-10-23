"""Workflow diagram for docs/ed-stats/PLAN.md (hand-authored SVG; render with rsvg-convert)."""
import os
D = os.path.dirname(os.path.abspath(__file__))
C = {"reused": ("#f1f0ec", "#8a8880"), "new": ("#e3f4ea", "#1f8a4c"), "mod": ("#efe6f7", "#7a4aa8")}
W, H = 980, 640
boxes = [  # x, y, w, h, kind, title, sub
    (30, 70, 200, 70, "reused", "read batch", "worker_pipeline step 0"),
    (270, 70, 200, 70, "mod", "worker_for_seq", "+ refmap_ed per read"),
    (510, 40, 200, 60, "reused", "placement", "refmap_query, pairing"),
    (750, 40, 200, 60, "reused", "stdout, PS4G, npy", "unchanged (byte-identical)"),
    (510, 130, 200, 70, "new", "tiles", "4 disjoint exact seeds\nrepeats: extend right"),
    (510, 230, 200, 80, "new", "best-first search", "match/sub/ins/del, 0..8\nprune: done rows, tiles,\nU, --ed-budget"),
    (510, 340, 200, 70, "new", "locate end states", "rb3_ssa_multi\nmin per assembly"),
    (270, 450, 200, 70, "mod", "write_refmap", "+ refmap_ed_write"),
    (30, 560, 200, 60, "new", "--ed-stats", "TSV or uint8 .npy"),
    (270, 560, 200, 60, "new", "--ed-hist", "best / asm / wrong / nbest"),
    (510, 560, 200, 60, "reused", "2-Poisson fit", "errprof.c (rb3_err_fit)"),
]
arrows = [((230, 105), (270, 105)), ((470, 95), (510, 70)), ((710, 70), (750, 70)), ((470, 115), (510, 165)),
          ((610, 200), (610, 230)), ((610, 310), (610, 340)), ((510, 375), (470, 485)), ((370, 520), (130, 560)),
          ((370, 520), (370, 560)), ((510, 590), (470, 590))]
o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Liberation Sans">',
     f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
     '<defs><marker id="a" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#5f5e5a"/></marker></defs>',
     '<text x="30" y="30" font-size="18" fill="#141413">refmap --ed-stats / --ed-hist (branch local/ed-stats)</text>']
for (x1, y1), (x2, y2) in arrows:
    o.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#5f5e5a" stroke-width="1.5" marker-end="url(#a)"/>')
for x, y, w, h, k, t, sub in boxes:
    f, s = C[k]
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{f}" stroke="{s}" stroke-width="2"/>')
    o.append(f'<text x="{x+12}" y="{y+22}" font-size="14" font-weight="bold" fill="#141413">{t}</text>')
    for i, line in enumerate(sub.split("\n")):
        o.append(f'<text x="{x+12}" y="{y+40+i*14}" font-size="11" fill="#3d3c38" font-family="Liberation Mono">{line}</text>')
lx = 750
for i, (k, lab) in enumerate((("reused", "reused, unchanged"), ("new", "new"), ("mod", "modified"))):
    f, s = C[k]
    o.append(f'<rect x="{lx}" y="{150+i*30}" width="22" height="16" rx="3" fill="{f}" stroke="{s}" stroke-width="2"/>')
    o.append(f'<text x="{lx+30}" y="{163+i*30}" font-size="12" fill="#141413">{lab}</text>')
o.append('<rect x="740" y="290" width="215" height="230" rx="8" fill="#fbfaf7" stroke="#c9c7bf"/>')
for i, line in enumerate(["Acceptance", "flags off: byte-identical to", "  ropebwt3_ep_a0fec31 (SE + PE)", "ED on: refmap outputs unchanged",
                          "unit + integration + valgrind", "agreement vs minimap2 (10k, 1.5M)", "overhead on Ia453 0.1x",
                          "A/B/C figure side by side"]):
    bold = 'font-weight="bold"' if i == 0 else ''
    o.append(f'<text x="752" y="{315+i*26}" font-size="{13 if i == 0 else 11}" {bold} fill="#141413">{line}</text>')
o.append("</svg>")
open(f"{D}/plan_diagram.svg", "w").write("\n".join(o))
