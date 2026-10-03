"""Shape (VSEPR), bond polarity, and molecular polarity, worked out from a Lewis structure.

Everything here is computed from the same SMILES the Lewis drawings use, so a
shape name or polarity call on the page can't disagree with the structure.
"""
import html
import re
import math

import lewis

# Electronegativity values as printed in most introductory textbooks (Pauling scale, one decimal).
EN = {"H": 2.1, "C": 2.5, "N": 3.0, "O": 3.5, "F": 4.0, "Na": 0.9, "P": 2.1, "S": 2.5,
      "Cl": 3.0, "Br": 2.8, "I": 2.5, "Xe": 2.6}
POLAR_MIN = 0.5   # ΔEN at or above this: polar covalent
IONIC_MIN = 2.0   # ΔEN at or above this: ionic

# (electron groups, lone pairs) -> (electron geometry, molecular geometry, bond angle)
GEOMETRY = {
    (2, 0): ("linear", "linear", "180°"),
    (3, 0): ("trigonal planar", "trigonal planar", "120°"),
    (3, 1): ("trigonal planar", "bent", "a little less than 120°"),
    (4, 0): ("tetrahedral", "tetrahedral", "109.5°"),
    (4, 1): ("tetrahedral", "trigonal pyramidal", "about 107°"),
    (4, 2): ("tetrahedral", "bent", "about 104.5°"),
    (5, 0): ("trigonal bipyramidal", "trigonal bipyramidal", "90° and 120°"),
    (5, 1): ("trigonal bipyramidal", "seesaw", "a little less than 90° and 120°"),
    (5, 2): ("trigonal bipyramidal", "T-shaped", "a little less than 90°"),
    (5, 3): ("trigonal bipyramidal", "linear", "180°"),
    (6, 0): ("octahedral", "octahedral", "90°"),
    (6, 1): ("octahedral", "square pyramidal", "a little less than 90°"),
    (6, 2): ("octahedral", "square planar", "90°"),
}


def bond_type(a, b):
    d = abs(EN[a] - EN[b])
    return "ionic" if d >= IONIC_MIN else "polar" if d >= POLAR_MIN else "nonpolar"


def center_of(m):
    """The central atom: the one bonded to the most atoms (ties go to the earlier atom)."""
    atoms = [a for a in m.GetAtoms() if a.GetAtomicNum()]
    return max(atoms, key=lambda a: (a.GetDegree(), -a.GetIdx())).GetIdx()


def vsepr(m, idx=None):
    idx = center_of(m) if idx is None else idx
    a = m.GetAtomWithIdx(idx)
    groups, lp = a.GetDegree() + lewis._lone_pairs(a), lewis._lone_pairs(a)
    if (groups, lp) not in GEOMETRY:
        raise ValueError(f"no VSEPR entry for {groups} groups, {lp} lone pairs")
    eg, mg, angle = GEOMETRY[(groups, lp)]
    return dict(center=idx, symbol=a.GetSymbol(), groups=groups, lone_pairs=lp,
                electron=eg, molecular=mg, angle=angle)


def partials(m):
    """δ+ / δ− for each atom in at least one polar bond. Returns {} for ions."""
    if any(a.GetFormalCharge() for a in m.GetAtoms()):
        return {}
    out = {}
    for a in m.GetAtoms():
        pull = 0.0
        for n in a.GetNeighbors():
            d = EN[a.GetSymbol()] - EN[n.GetSymbol()]
            if abs(d) >= POLAR_MIN:
                pull += d
        if pull:
            out[a.GetIdx()] = "δ−" if pull > 0 else "δ+"
    return out


def polarity(m):
    """('polar' | 'nonpolar', reason), using the rules taught in this course."""
    polar_bonds = [b for b in m.GetBonds()
                   if bond_type(b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()) != "nonpolar"]
    if not polar_bonds:
        return "nonpolar", "no polar bonds"
    centers = [a for a in m.GetAtoms() if a.GetDegree() > 1]
    if not centers:
        return "polar", "one polar bond"
    if len(centers) > 1:
        return "polar", "polar bonds spread over more than one central atom"
    c = centers[0]
    outer = {(n.GetSymbol(), int(m.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx()).GetBondTypeAsDouble()))
             for n in c.GetNeighbors()}
    if lewis._lone_pairs(c) == 0 and len(outer) == 1:
        return "nonpolar", "identical polar bonds arranged symmetrically, so they cancel"
    if lewis._lone_pairs(c):
        return "polar", "lone pairs on the central atom make the shape lopsided"
    return "polar", "different atoms around the center, so the bond polarities don't cancel"


# ------------------------------------------------------------------ 3D shape drawings
# Hand-placed 2D projections: (kind, x, y) with kind P (in the page), W (wedge, toward
# you), D (dash, away). y is up. One list per (groups, lone pairs); lone pairs take the
# slots listed in LP_SLOTS, atoms fill the rest in order.
SLOTS = {
    (2, 0): [("P", -1, 0), ("P", 1, 0)],
    (3, 0): [("P", 0, 1), ("P", -0.87, -0.5), ("P", 0.87, -0.5)],
    (3, 1): [("P", 0, 1), ("P", -0.87, -0.5), ("P", 0.87, -0.5)],
    (4, 0): [("P", 0, 1), ("P", -0.94, -0.34), ("W", 0.34, -0.94), ("D", 0.98, -0.17)],
    (4, 1): [("P", 0, 1), ("P", -0.94, -0.34), ("W", 0.34, -0.94), ("D", 0.98, -0.17)],
    (4, 2): [("W", 0.62, 0.78), ("D", -0.62, 0.78), ("P", -0.82, -0.57), ("P", 0.82, -0.57)],
    (6, 0): [("P", 0, 1), ("P", 0, -1), ("W", -0.85, -0.45), ("W", 0.85, -0.45), ("D", 0.85, 0.45), ("D", -0.85, 0.45)],
    (6, 1): [("P", 0, 1), ("P", 0, -1), ("W", -0.85, -0.45), ("W", 0.85, -0.45), ("D", 0.85, 0.45), ("D", -0.85, 0.45)],
}
LP_SLOTS = {(3, 1): [0], (4, 1): [0], (4, 2): [0, 1], (6, 1): [1]}
L = 50  # bond length in pixels


def shape_svg(uid, m, title, idx=None):
    g = vsepr(m, idx)
    key = (g["groups"], g["lone_pairs"])
    if key not in SLOTS:
        raise ValueError(f"{title}: no shape drawing for {key}")
    slots = SLOTS[key]
    lp_slots = LP_SLOTS.get(key, [])
    free = [s for k, s in enumerate(slots) if k not in lp_slots]
    c = m.GetAtomWithIdx(g["center"])
    nbrs = sorted(c.GetNeighbors(), key=lambda n: (
        -m.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx()).GetBondTypeAsDouble(), n.GetSymbol() == "H", n.GetIdx()))
    parts, labels = [], []
    pts = [(0, 0)]

    def lab(sym, x, y):
        labels.append(f'<text x="{x:.1f}" y="{y:.1f}" dy="0.35em" style="fill:var({lewis.COLOR.get(sym, "--mol-c")})">{sym}</text>')

    def charge(atom, x, y, used, col):
        fc = atom.GetFormalCharge()
        if not fc:
            return
        tag = lewis.charge_text(fc)
        cand = [k * 15 for k in range(24)]
        best = max(cand, key=lambda a: (min((lewis._adiff(a, u) for u in used), default=180), -lewis._adiff(a, -45)))
        t = math.radians(best)
        rc = lewis._wide(atom.GetSymbol(), t, lewis.R_CHARGE + (4 if len(tag) > 1 else 0))
        cx, cy = x + rc * math.cos(t), y + rc * math.sin(t)
        labels.append(f'<text class="fc" x="{cx:.1f}" y="{cy:.1f}" dy="0.35em" style="fill:{col}">{tag}</text>')
        pts.append((cx, cy))

    lab(c.GetSymbol(), 0, 0)
    for n, (kind, sx, sy) in zip(nbrs, free):
        x, y = sx * L, -sy * L
        pts.append((x, y))
        d = math.hypot(x, y)
        ux, uy = x / d, y / d
        r0, r1 = lewis._wide(c.GetSymbol(), math.atan2(uy, ux), 11), lewis._wide(n.GetSymbol(), math.atan2(uy, ux), 11)
        x0, y0, x1, y1 = ux * r0, uy * r0, x - ux * r1, y - uy * r1
        px, py = -uy, ux
        if kind == "P":
            order = int(m.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx()).GetBondTypeAsDouble())
            for o in {1: [0], 2: [-3.2, 3.2], 3: [-5, 0, 5]}[order]:
                parts.append(f'<line x1="{x0 + px * o:.1f}" y1="{y0 + py * o:.1f}" x2="{x1 + px * o:.1f}" y2="{y1 + py * o:.1f}"/>')
        elif kind == "W":
            w = 4.2
            parts.append(f'<path class="wedge" d="M{x0:.1f},{y0:.1f} L{x1 + px * w:.1f},{y1 + py * w:.1f} '
                         f'L{x1 - px * w:.1f},{y1 - py * w:.1f} Z"/>')
        else:
            for k in range(1, 7):
                t = k / 7
                hx, hy = x0 + (x1 - x0) * t, y0 + (y1 - y0) * t
                w = 0.8 + 4 * t
                parts.append(f'<line class="hash" x1="{hx + px * w:.1f}" y1="{hy + py * w:.1f}" '
                             f'x2="{hx - px * w:.1f}" y2="{hy - py * w:.1f}"/>')
        lab(n.GetSymbol(), x, y)
    # Outer atoms: draw their lone pairs and formal charges the same way the Lewis drawings do,
    # so every C, N, O, and F still shows its full octet in the 3D picture.
    for n, (kind, sx, sy) in zip(nbrs, free):
        x, y = sx * L, -sy * L
        toward = math.degrees(math.atan2(-y, -x))
        ncolor = f"var({lewis.COLOR.get(n.GetSymbol(), '--mol-c')})"
        slots_n = lewis.spread([toward], lewis._lone_pairs(n), prefer=toward + 180)
        for t in slots_n:
            parts.append(lewis._dot_pair(x, y, t, ncolor, n.GetSymbol()))
            pts.append((x + 20 * math.cos(math.radians(t)), y + 20 * math.sin(math.radians(t))))
        charge(n, x, y, [toward] + slots_n, ncolor)
    color = f"var({lewis.COLOR.get(c.GetSymbol(), '--mol-c')})"
    for k in lp_slots:
        kind, sx, sy = slots[k]
        x, y = sx * L * 0.55, -sy * L * 0.55
        d = math.hypot(x, y)
        px, py = -y / d * 3.6, x / d * 3.6
        parts.append(f'<circle cx="{x + px:.1f}" cy="{y + py:.1f}" r="2.4" fill="{color}"/>'
                     f'<circle cx="{x - px:.1f}" cy="{y - py:.1f}" r="2.4" fill="{color}"/>')
        pts.append((x, y))
    center_used = [math.degrees(math.atan2(-sy, sx)) for _, sx, sy in slots]
    charge(c, 0, 0, center_used, color)
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    x0, x1, y0, y1 = min(xs) - 22, max(xs) + 22, min(ys) - 20, max(ys) + 20
    w, h = x1 - x0, y1 - y0
    desc = (f"Three-dimensional shape of {title}: {g['molecular']}. The central {lewis.WORD[g['symbol']]} has "
            f"{g['groups']} electron groups ({g['lone_pairs']} lone pair{'s' if g['lone_pairs'] != 1 else ''}). "
            "Each outer atom shows its own lone pairs. Solid wedges point toward you; dashed bonds point away.")
    return (f'<svg class="lewis shape" viewBox="{x0:.1f} {y0:.1f} {w:.1f} {h:.1f}" width="{w * lewis.DISPLAY:.0f}" '
            f'height="{h * lewis.DISPLAY:.0f}" role="img" aria-labelledby="t-{uid}" xmlns="http://www.w3.org/2000/svg">'
            f'<title id="t-{uid}">{html.escape(desc)}</title><g class="bonds">{"".join(parts)}</g>'
            f'<g class="atoms">{"".join(labels)}</g></svg>'), g


# ------------------------------------------------------------------ intermolecular force diagrams
def _water(cx, cy, toward, h_in, parts, labels, atoms):
    """A water molecule centered on its O at (cx, cy).

    toward: angle (degrees, screen coordinates) pointing at the neighbor it interacts with.
    h_in:   True if an H points at the neighbor, False if the O does.
    atoms:  list that collects (x, y, symbol) so every atom can get its partial charge later.
    """
    half = 52.25  # half of 104.5°
    hs = [toward, toward + 2 * half] if h_in else [toward + 180 - half, toward + 180 + half]
    r = 30
    for t in hs:
        a = math.radians(t)
        x, y = cx + r * math.cos(a), cy + r * math.sin(a)
        ux, uy = math.cos(a), math.sin(a)
        parts.append(f'<line x1="{cx + ux * 10:.1f}" y1="{cy + uy * 10:.1f}" x2="{x - ux * 8:.1f}" y2="{y - uy * 8:.1f}"/>')
        labels.append(f'<text x="{x:.1f}" y="{y:.1f}" dy="0.35em" style="fill:var(--mol-c)">H</text>')
        atoms.append((x, y, "H"))
    labels.append(f'<text x="{cx:.1f}" y="{cy:.1f}" dy="0.35em" style="fill:var(--mol-o)">O</text>')
    atoms.append((cx, cy, "O"))
    return [(cx + r * math.cos(math.radians(t)), cy + r * math.sin(math.radians(t))) for t in hs]


_LINE = re.compile(r'<line[^>]*x1="([-\d.]+)" y1="([-\d.]+)" x2="([-\d.]+)" y2="([-\d.]+)"')


def _partial_charges(parts, labels, atoms):
    """Label every O with δ− and every H with δ+, in the most open direction around each atom
    (away from its bonds and from the dotted attraction lines)."""
    segs = [tuple(map(float, m.groups())) for m in map(_LINE.search, parts) if m]
    for x, y, sym in atoms:
        occupied = []
        for x1, y1, x2, y2 in segs:
            for (ex, ey), (fx, fy) in (((x1, y1), (x2, y2)), ((x2, y2), (x1, y1))):
                if math.hypot(ex - x, ey - y) < 16:
                    occupied.append(math.degrees(math.atan2(fy - y, fx - x)))
        t = math.radians(lewis.spread(occupied, 1, prefer=-90)[0])
        rr = 19 if sym == "O" else 17
        tag, col = ("δ−", "--mol-o") if sym == "O" else ("δ+", "--mol-c")
        labels.append(f'<text class="pc" x="{x + rr * math.cos(t):.1f}" y="{y + rr * math.sin(t):.1f}" dy="0.35em" '
                      f'style="fill:var({col})">{tag}</text>')


def _svg(uid, parts, labels, box, desc):
    x0, y0, x1, y1 = box
    w, h = x1 - x0, y1 - y0
    return (f'<svg class="lewis imf" viewBox="{x0} {y0} {w} {h}" width="{w * lewis.DISPLAY:.0f}" '
            f'height="{h * lewis.DISPLAY:.0f}" role="img" aria-labelledby="t-{uid}" xmlns="http://www.w3.org/2000/svg">'
            f'<title id="t-{uid}">{html.escape(desc)}</title><g class="bonds">{"".join(parts)}</g>'
            f'<g class="atoms">{"".join(labels)}</g></svg>')


def hydration_svg():
    """Na+ with water O atoms facing it; Cl− with water H atoms facing it."""
    parts, labels, atoms = [], [], []

    def dotted(x0, y0, x1, y1):
        parts.append(f'<line class="imfline" x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}"/>')

    d = 72
    for k, ang in enumerate([0, 90, 180, 270]):
        a = math.radians(ang)
        cx, cy = d * math.cos(a), d * math.sin(a)
        _water(cx, cy, ang + 180, False, parts, labels, atoms)
        dotted(cx - 11 * math.cos(a), cy - 11 * math.sin(a), 20 * math.cos(a), 20 * math.sin(a))
    labels.append('<text class="ionlabel" x="0" y="0" dy="0.35em" style="fill:var(--mol-c)">Na⁺</text>')
    ox = 280
    for k, ang in enumerate([45, 135, 225, 315]):
        a = math.radians(ang)
        cx, cy = ox + 80 * math.cos(a), 80 * math.sin(a)
        _water(cx, cy, ang + 180, True, parts, labels, atoms)
        hx, hy = ox + 50 * math.cos(a), 50 * math.sin(a)
        dotted(hx - 9 * math.cos(a), hy - 9 * math.sin(a), ox + 19 * math.cos(a), 19 * math.sin(a))
    labels.append(f'<text class="ionlabel" x="{ox}" y="0" dy="0.35em" style="fill:var(--mol-x)">Cl⁻</text>')
    _partial_charges(parts, labels, atoms)
    desc = ("Ion–dipole attractions. Left: a sodium ion surrounded by four water molecules, each with its "
            "partly negative oxygen end facing the positive ion. Right: a chloride ion surrounded by four water "
            "molecules, each with one partly positive hydrogen pointing at the negative ion. Dotted lines show "
            "the attractions. Every oxygen is labeled partial negative and every hydrogen partial positive.")
    return _svg("hydration", parts, labels, (-125, -125, 400, 125), desc)


def hbond_svg():
    """Three water molecules joined by hydrogen bonds (dotted)."""
    parts, labels, atoms = [], [], []
    centers = [(0, 0), (104, 0), (52, 92)]
    # molecule 0 donates an H to molecule 1; molecule 1 donates an H to molecule 2
    hs0 = _water(*centers[0], 0, True, parts, labels, atoms)
    hs1 = _water(*centers[1], math.degrees(math.atan2(92, 52 - 104)), True, parts, labels, atoms)
    _water(*centers[2], -90, False, parts, labels, atoms)
    for (hx, hy), (ox, oy) in [(hs0[0], centers[1]), (hs1[0], centers[2])]:
        dx, dy = ox - hx, oy - hy
        dd = math.hypot(dx, dy)
        parts.append(f'<line class="imfline" x1="{hx + dx / dd * 9:.1f}" y1="{hy + dy / dd * 9:.1f}" '
                     f'x2="{ox - dx / dd * 11:.1f}" y2="{oy - dy / dd * 11:.1f}"/>')
    _partial_charges(parts, labels, atoms)
    desc = ("Hydrogen bonds in water. Three water molecules: a hydrogen on one molecule points at the oxygen of the "
            "next, joined by a dotted line representing the hydrogen bond, which is an attraction between molecules, "
            "not a covalent bond. Every oxygen is labeled partial negative and every hydrogen partial positive.")
    return _svg("hbond", parts, labels, (-45, -45, 160, 140), desc)
