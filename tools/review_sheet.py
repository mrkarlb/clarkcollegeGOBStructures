"""Build _site/structure_review.html: every structure on one page, with its
electron count around each atom, for chemistry review before publishing.

    python tools/review_sheet.py
"""
import html
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import build  # noqa: E402
import lewis  # noqa: E402

build.verify()
cells = []
for s in build.S:
    m = build.MOLS[s["id"]]
    atoms = []
    for a in m.GetAtoms():
        if not a.GetAtomicNum() or a.GetSymbol() == "H":
            continue
        e = 2 * lewis._bond_sum(a) + 2 * lewis._lone_pairs(a)
        fc = a.GetFormalCharge()
        atoms.append(f"{a.GetSymbol()}{a.GetIdx()}: {e} e⁻, {lewis._lone_pairs(a)} LP"
                     + (f", {lewis.charge_text(fc)}" if fc else ""))
    total = lewis.electron_count_rows(s["formula"])[2] if s["formula"] else "—"
    cells.append(
        f'<div class="cell"><div class="pic">{build.drawing(s["id"])}</div>'
        f'<b>{html.escape(s["id"])}</b> · {lewis.formula_html(s["formula"]) if s["formula"] else ""} '
        f'{html.escape(s["title"])}<br><small>{html.escape(s["smiles"])} · role {s["role"]} · '
        f'{total} valence e⁻<br>{"<br>".join(atoms)}'
        f'{"<br>OPSIN: " + html.escape(s["opsin"]) if s["opsin"] else ""}</small></div>')
tpl = open(os.path.join(build.ROOT, "template.html"), encoding="utf-8").read()
style = tpl[tpl.index("<style>"):tpl.index("</style>") + 8]
page = (f'<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><title>Structure review</title>{style}'
        '<style>.grid{display:flex;flex-wrap:wrap;gap:12px;padding:16px}.cell{background:var(--surface);'
        'border:1px solid var(--border);border-radius:6px;padding:8px;width:260px;font-size:13px}</style></head>'
        f'<body><h1 style="padding:0 16px">Structure review: {len(cells)} structures</h1>'
        f'<div class="grid">{"".join(cells)}</div></body></html>')
os.makedirs(build.OUT, exist_ok=True)
open(os.path.join(build.OUT, "structure_review.html"), "w", encoding="utf-8").write(page)
print("Wrote _site/structure_review.html")
