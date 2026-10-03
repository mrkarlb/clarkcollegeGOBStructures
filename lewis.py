"""Lewis structures: electron bookkeeping and theme-aware SVG drawings.

Every structure is entered as SMILES with formal charges ([NH4+], [C-]#[O+]).
From that alone this module works out each atom's lone pairs, checks the
electron count around every atom, and draws the structure with all atoms,
bonds, lone pairs, and formal charges shown.

Colors are CSS custom properties (--mol-c, --mol-o, ...), defined for light
and dark themes in template.html.
"""
import html
import math
import re

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point2D

rdDepictor.SetPreferCoordGen(True)

# Valence electrons (group number for main-group elements)
VALENCE = {"H": 1, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Na": 1, "P": 5, "S": 6,
           "Cl": 7, "Br": 7, "I": 7, "Xe": 8}
OCTET_ONLY = {"C", "N", "O", "F"}            # second shell: s and p only, never more than 8
CAN_EXPAND = {"P", "S", "Cl", "Br", "I", "Xe"}  # larger atoms: room for more than 8
COLOR = {"C": "--mol-c", "H": "--mol-c", "O": "--mol-o", "N": "--mol-n", "S": "--mol-s",
         "F": "--mol-x", "Cl": "--mol-x", "Br": "--mol-x", "I": "--mol-x", "P": "--mol-p",
         "Xe": "--mol-xe", "B": "--mol-c", "Na": "--mol-c"}
WORD = {"H": "hydrogen", "C": "carbon", "N": "nitrogen", "O": "oxygen", "F": "fluorine",
        "P": "phosphorus", "S": "sulfur", "Cl": "chlorine", "Br": "bromine", "I": "iodine",
        "Xe": "xenon", "B": "boron", "Na": "sodium"}
BOND_WORD = {1: "single", 2: "double", 3: "triple"}

SCALE = 50       # pixels per layout unit (one bond length)
R_LABEL = 11     # bonds stop this far from a labeled atom's center
R_DOT = 15       # distance from atom center to its lone-pair dots
DOT_SEP = 3.6    # half the spacing between the two dots of a pair
DOT_R = 2.3
R_CHARGE = 24
DISPLAY = 1.15  # on-page size relative to the drawing


# ------------------------------------------------------------------ formulas
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
SUP = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")


def parse_formula(f):
    """'CO3^2-' -> ({'C': 1, 'O': 3}, -2)"""
    body, _, chg = f.partition("^")
    counts = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", body):
        counts[el] = counts.get(el, 0) + (int(n) if n else 1)
    charge = 0
    if chg:
        m = re.fullmatch(r"(\d*)([+-])", chg)
        if not m:
            raise ValueError(f"bad charge in formula {f}")
        charge = (int(m.group(1)) if m.group(1) else 1) * (1 if m.group(2) == "+" else -1)
    return counts, charge


def formula_html(f):
    """'CO3^2-' -> 'CO₃²⁻' (plain Unicode so it reads correctly everywhere)"""
    body, _, chg = f.partition("^")
    out = re.sub(r"(\d+)", lambda m: m.group(1).translate(SUB), body)
    return out + chg.replace("-", "−").translate(SUP).replace("−", "⁻") if chg else out


def charge_text(q):
    if q == 0:
        return ""
    n = "" if abs(q) == 1 else str(abs(q))
    return n + ("+" if q > 0 else "−")


# ------------------------------------------------------------------ molecule
def make_mol(smiles, coords=None, like=None, rotate=0, flip=False):
    """Parse SMILES, add hydrogens, and lay the molecule out in 2D.

    coords: {atom index: (x, y)} to pin heavy atoms (y up, one unit per bond)
    like:   another laid-out molecule with the same atoms in the same order
            (resonance forms share one layout)
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"RDKit could not read SMILES {smiles}")
    # RDKit quietly rewrites some structures (a five-bond N becomes N+ / O-).
    # Refuse those, so what's drawn is exactly what was written.
    raw = Chem.MolFromSmiles(smiles, sanitize=False)
    as_written = [(a.GetFormalCharge(), sorted(b.GetBondType() for b in a.GetBonds())) for a in raw.GetAtoms()]
    as_read = [(a.GetFormalCharge(), sorted(b.GetBondType() for b in a.GetBonds())) for a in m.GetAtoms()]
    if as_written != as_read:
        raise ValueError(f"RDKit had to rewrite {smiles} as {Chem.MolToSmiles(m)}; "
                         "check the bonds and charges")
    m = Chem.AddHs(m)
    if like is not None:
        if [a.GetSymbol() for a in m.GetAtoms()] != [a.GetSymbol() for a in like.GetAtoms()]:
            raise ValueError(f"{smiles}: atoms do not match the structure it borrows its layout from")
        conf = Chem.Conformer(m.GetNumAtoms())
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, like.GetConformer().GetAtomPosition(i))
        m.AddConformer(conf, assignId=True)
        return m
    cmap = {i: Point2D(x, y) for i, (x, y) in (coords or {}).items()}
    rdDepictor.Compute2DCoords(m, coordMap=cmap) if cmap else rdDepictor.Compute2DCoords(m)
    conf = m.GetConformer()
    if coords:
        # keep pinned atoms exactly where they were asked to be
        for i, (x, y) in coords.items():
            conf.SetAtomPosition(i, (x, y, 0))
        _place_hydrogens(m, set(coords))
    a = math.radians(rotate)
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        x, y = p.x * math.cos(a) - p.y * math.sin(a), p.x * math.sin(a) + p.y * math.cos(a)
        conf.SetAtomPosition(i, (-x if flip else x, y, 0))
    return m


def _place_hydrogens(m, pinned):
    """Put each unpinned H one bond length out, in the widest open direction around its atom."""
    conf = m.GetConformer()
    for atom in m.GetAtoms():
        hs = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "H" and n.GetIdx() not in pinned]
        if not hs:
            continue
        c = conf.GetAtomPosition(atom.GetIdx())
        fixed = [math.degrees(math.atan2(conf.GetAtomPosition(n.GetIdx()).y - c.y,
                                         conf.GetAtomPosition(n.GetIdx()).x - c.x))
                 for n in atom.GetNeighbors() if n.GetIdx() not in hs]
        lone = _lone_pairs(atom)
        angles = spread(fixed, len(hs) + lone, prefer=-90)
        # hydrogens take the slots farthest from the existing bonds; lone pairs get the rest
        angles.sort(key=lambda t: -min((_adiff(t, f) for f in fixed), default=0))
        for h, t in zip(hs, angles[:len(hs)]):
            conf.SetAtomPosition(h,(c.x + math.cos(math.radians(t)), c.y + math.sin(math.radians(t)), 0))


# ------------------------------------------------------------------ bookkeeping
def _bond_sum(atom):
    return sum(int(b.GetBondTypeAsDouble()) for b in atom.GetBonds())


def _lone_pairs(atom):
    v = VALENCE.get(atom.GetSymbol())
    if v is None:
        return 0
    return (v - atom.GetFormalCharge() - _bond_sum(atom)) // 2


def audit(m, formula, label):
    """Return a list of problems with this structure (empty list = it checks out)."""
    problems = []
    counts, charge = parse_formula(formula)
    have = {}
    for a in m.GetAtoms():
        if a.GetAtomicNum() == 0:
            continue
        have[a.GetSymbol()] = have.get(a.GetSymbol(), 0) + 1
    if have != counts:
        problems.append(f"{label}: SMILES gives {have}, formula {formula} says {counts}")
    total_fc = sum(a.GetFormalCharge() for a in m.GetAtoms())
    if total_fc != charge:
        problems.append(f"{label}: formal charges add to {total_fc}, but the charge is {charge}")
    for a in m.GetAtoms():
        sym = a.GetSymbol()
        if a.GetAtomicNum() == 0:
            continue
        if sym not in VALENCE:
            problems.append(f"{label}: no valence-electron count for {sym}")
            continue
        free = VALENCE[sym] - a.GetFormalCharge() - _bond_sum(a)
        where = f"{label}: {sym} (atom {a.GetIdx()})"
        if free < 0:
            problems.append(f"{where} has more bonds than its electrons and charge allow")
            continue
        if free % 2:
            problems.append(f"{where} has an unpaired electron")
            continue
        shell = 2 * _bond_sum(a) + free
        if sym == "H":
            if shell != 2:
                problems.append(f"{where} has {shell} electrons; H needs 2")
        elif sym in OCTET_ONLY:
            if shell != 8:
                problems.append(f"{where} has {shell} electrons; {sym} must have exactly 8")
        elif sym in CAN_EXPAND:
            groups = a.GetDegree() + free // 2
            if shell < 8:
                problems.append(f"{where} has {shell} electrons; expected at least 8")
            if groups > 6:
                problems.append(f"{where} has {groups} electron groups; six is the most in this course")
        elif shell != 8:
            problems.append(f"{where} has {shell} electrons; expected 8")
    # the whole-molecule count, done the way students do it
    expected = sum(VALENCE[el] * n for el, n in counts.items() if el in VALENCE) - charge
    bonds =sum(int(b.GetBondTypeAsDouble()) for b in m.GetBonds() if 0 not in
                (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()))
    drawn = 2 * bonds + sum(2 * _lone_pairs(a) for a in m.GetAtoms() if a.GetAtomicNum())
    if drawn != expected:
        problems.append(f"{label}: {drawn} electrons drawn, but the formula gives {expected}")
    return problems


def expanded_atoms(m):
    return [a.GetSymbol() for a in m.GetAtoms()
            if a.GetAtomicNum() and 2 * _bond_sum(a) + 2 * _lone_pairs(a) > 8]


def electron_count_rows(formula):
    counts, charge = parse_formula(formula)
    rows = [(el, n, VALENCE[el], n * VALENCE[el]) for el, n in counts.items()]
    total = sum(r[3] for r in rows) - charge
    return rows, charge, total


# ------------------------------------------------------------------ geometry helpers
def _adiff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def spread(occupied, k, prefer=-90):
    """Choose k directions (degrees) that sit in the open space around an atom.

    Free slots go into the widest gaps between the occupied directions, shared
    out in proportion to gap size, and evenly spaced within each gap.
    """
    if k <= 0:
        return []
    if not occupied:
        return [(prefer + i * 360 / k) % 360 for i in range(k)] if k != 4 \
            else [prefer % 360, (prefer + 90) % 360, (prefer + 180) % 360, (prefer + 270) % 360]
    occ = sorted(o % 360 for o in occupied)
    gaps = [(occ[i], (occ[(i + 1) % len(occ)] - occ[i]) % 360 or 360) for i in range(len(occ))]
    share = [0] * len(gaps)
    for _ in range(k):
        i = max(range(len(gaps)), key=lambda j: (gaps[j][1] / (share[j] + 1), -j))
        share[i] += 1
    out = []
    for (start, size), n in zip(gaps, share):
        out += [(start + size * (j + 1) / (n + 1)) % 360 for j in range(n)]
    return out


# ------------------------------------------------------------------ drawing
def _var(sym):
    return f"var({COLOR.get(sym, '--mol-c')})"


def _wide(sym, t, base):
    """Distance from an atom's center that clears its label (two-letter symbols are wider)."""
    return base + (5 if len(sym) > 1 else 0) * abs(math.cos(t))


def _dot_pair(x, y, ang, color, sym="C"):
    t = math.radians(ang)
    r = _wide(sym, t, R_DOT)
    cx, cy = x + r * math.cos(t), y + r * math.sin(t)
    px, py = -math.sin(t) * DOT_SEP, math.cos(t) * DOT_SEP
    return (f'<circle cx="{cx + px:.1f}" cy="{cy + py:.1f}" r="{DOT_R}" fill="{color}"/>'
            f'<circle cx="{cx - px:.1f}" cy="{cy - py:.1f}" r="{DOT_R}" fill="{color}"/>')


def _dot_single(x, y, ang, color, sym="C"):
    t = math.radians(ang)
    r = _wide(sym, t, R_DOT)
    return f'<circle cx="{x + r * math.cos(t):.1f}" cy="{y + r * math.sin(t):.1f}" r="{DOT_R}" fill="{color}"/>'


def _sigma_dots(atom):
    """Electrons left on an atom after one single bond to each neighbor,
    starting from its Lewis symbol: (pairs, single electrons)."""
    v = VALENCE[atom.GetSymbol()]
    singles = v if v <= 4 else 8 - v
    pairs = (v - singles) // 2
    singles -= atom.GetDegree()
    if singles < 0:
        raise ValueError(f"{atom.GetSymbol()}: more bonds than single electrons; "
                         "step drawings only work for molecules that need no lone pair to be shared")
    return pairs, singles


def lewis_svg(uid, m, alt, stage="final", brackets=True, notes=None):
    """Draw a laid-out molecule.

    stage: 'skeleton' (single bonds only), 'dots' (single bonds plus every
    atom's remaining electrons, as in Step 2), or 'final'.
    """
    conf = m.GetConformer()
    pos = {a.GetIdx(): (conf.GetAtomPosition(a.GetIdx()).x * SCALE, -conf.GetAtomPosition(a.GetIdx()).y * SCALE)
           for a in m.GetAtoms()}
    labeled = {a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()}
    bond_parts, dot_parts, text_parts = [], [], []
    marks = []  # every drawn point, for the bounding box

    def ang(i, j):
        (x1, y1), (x2, y2) = pos[i], pos[j]
        return math.degrees(math.atan2(y2 - y1, x2 - x1))

    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        (x1, y1), (x2, y2) = pos[i], pos[j]
        d = math.hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / d, (y2 - y1) / d
        t = math.atan2(uy, ux)
        r1 = _wide(m.GetAtomWithIdx(i).GetSymbol(), t, R_LABEL) if i in labeled else 0
        r2 = _wide(m.GetAtomWithIdx(j).GetSymbol(), t, R_LABEL) if j in labeled else 0
        order = 1 if stage != "final" else int(b.GetBondTypeAsDouble())
        offs = {1: [0], 2: [-3.2, 3.2], 3: [-5, 0, 5]}[order]
        for o in offs:
            px, py = -uy * o, ux * o
            bond_parts.append(f'<line x1="{x1 + ux * r1 + px:.1f}" y1="{y1 + uy * r1 + py:.1f}" '
                              f'x2="{x2 - ux * r2 + px:.1f}" y2="{y2 - uy * r2 + py:.1f}"/>')

    for a in m.GetAtoms():
        i, sym = a.GetIdx(), a.GetSymbol()
        x, y = pos[i]
        if not a.GetAtomicNum():
            marks.append((x, y))
            continue
        color = _var(sym)
        text_parts.append(f'<text x="{x:.1f}" y="{y:.1f}" dy="0.35em" style="fill:{color}">{sym}</text>')
        hw = 15 if len(sym) > 1 else 9
        marks += [(x - hw, y - 10), (x + hw, y + 10)]
        bonds = [ang(i, n.GetIdx()) for n in a.GetNeighbors()]
        if stage == "skeleton":
            continue
        if stage == "dots":
            pairs, singles = _sigma_dots(a)
            # a single electron that will become part of a double or triple bond
            # sits right beside that bond, so the pairing in Step 3 is easy to see
            beside = []
            for b in a.GetBonds():
                extra = int(b.GetBondTypeAsDouble()) - 1
                other = b.GetOtherAtomIdx(i)
                base = ang(i, other)
                side = 1 if i < other else -1  # both atoms' electrons land on the same side of the bond
                beside += [base + 52 * side * (1 if k % 2 == 0 else -1) for k in range(extra)]
            beside = beside[:singles]
            free = spread(bonds + beside, pairs + singles - len(beside))
            for t in beside + free[:singles - len(beside)]:
                dot_parts.append(_dot_single(x, y, t, color, sym))
            for t in free[singles - len(beside):]:
                dot_parts.append(_dot_pair(x, y, t, color, sym))
            used = bonds + beside + free
        else:
            lp = _lone_pairs(a)
            slots = spread(bonds, lp)
            for t in slots:
                dot_parts.append(_dot_pair(x, y, t, color, sym))
            used = bonds + slots
            fc = a.GetFormalCharge()
            tag = charge_text(fc) if fc else (notes or {}).get(i)
            if tag:
                cand = [c * 15 for c in range(24)]
                best = max(cand, key=lambda c: (min((_adiff(c, u) for u in used), default=180),
                                                -_adiff(c, -45)))
                t = math.radians(best)
                rc = _wide(sym, t, (17 if sym == "H" else R_CHARGE) + (4 if len(tag) > 1 else 0))
                cx, cy = x + rc * math.cos(t), y + rc * math.sin(t)
                cls = "fc" if fc else "pc"
                text_parts.append(f'<text class="{cls}" x="{cx:.1f}" y="{cy:.1f}" dy="0.35em" '
                                  f'style="fill:{color}">{tag}</text>')
                marks += [(cx - 9, cy - 8), (cx + 9, cy + 8)]
        for t in used:
            r = math.radians(t)
            rr = _wide(sym, r, R_DOT + 3)
            marks.append((x + rr * math.cos(r), y + rr * math.sin(r)))

    xs, ys = [p[0] for p in marks], [p[1] for p in marks]
    pad = 6
    x0, x1, y0, y1 = min(xs) - pad, max(xs) + pad, min(ys) - pad, max(ys) + pad
    net = sum(a.GetFormalCharge() for a in m.GetAtoms())
    extra = ""
    if stage == "final" and net and brackets:
        bx0, bx1 = x0 - 2, x1 + 2
        extra = (f'<path class="br" d="M{bx0 + 6:.1f},{y0:.1f} h-6 v{y1 - y0:.1f} h6 '
                 f'M{bx1 - 6:.1f},{y0:.1f} h6 v{y1 - y0:.1f} h-6"/>'
                 f'<text class="ion" x="{bx1 + 3:.1f}" y="{y0 + 8:.1f}">{charge_text(net)}</text>')
        x0, x1 = bx0 - 4, bx1 + 22
    w, h = x1 - x0, y1 - y0
    tid = f"t-{uid}"
    return (f'<svg class="lewis" viewBox="{x0:.1f} {y0:.1f} {w:.1f} {h:.1f}" width="{w * DISPLAY:.0f}" height="{h * DISPLAY:.0f}" '
            f'role="img" aria-labelledby="{tid}" xmlns="http://www.w3.org/2000/svg">'
            f'<title id="{tid}">{html.escape(alt)}</title>'
            f'<g class="bonds">{"".join(bond_parts)}</g><g>{"".join(dot_parts)}</g>'
            f'<g class="atoms">{"".join(text_parts)}</g>{extra}</svg>')


def describe(m, title, stage="final", notes=None):
    """Plain-language description for screen readers."""
    heavy = [a for a in m.GetAtoms() if a.GetAtomicNum() and a.GetSymbol() != "H"]
    if not heavy:
        heavy = [a for a in m.GetAtoms() if a.GetAtomicNum()]
    lead = {"skeleton": "First step for", "dots": "Second step for", "final": "Lewis structure of"}[stage]
    parts = [f"{lead} {title}."]
    for a in heavy:
        nbrs = {}
        for b in a.GetBonds():
            o = b.GetOtherAtom(a)
            if not o.GetAtomicNum():
                key = "an unnamed group"
            else:
                key = WORD[o.GetSymbol()]
            order = 1 if stage != "final" else int(b.GetBondTypeAsDouble())
            nbrs.setdefault((key, order), 0)
            nbrs[(key, order)] += 1
        bits = []
        for (key, order), n in nbrs.items():
            if key == "an unnamed group":
                bits.append(f"{n} {BOND_WORD[order]} bond{'s' if n > 1 else ''} to other atoms")
            else:
                bits.append(f"{n} {key}{'s' if n > 1 else ''} by {BOND_WORD[order]} bond{'s' if n > 1 else ''}")
        s = f"{WORD[a.GetSymbol()].capitalize()} is bonded to " + ", ".join(bits) if bits else WORD[a.GetSymbol()].capitalize()
        if stage == "final":
            lp = _lone_pairs(a)
            if lp:
                s += f", with {lp} lone pair{'s' if lp > 1 else ''}"
            fc = a.GetFormalCharge()
            if fc:
                s += f", and a formal charge of {'+' if fc > 0 else '−'}{abs(fc)}"
            elif notes and a.GetIdx() in notes:
                s += f", and a partial {'negative' if '−' in notes[a.GetIdx()] else 'positive'} charge"
        elif stage == "dots":
            p, sg = _sigma_dots(a)
            s += f", with {p} pair{'s' if p != 1 else ''} and {sg} single electron{'s' if sg != 1 else ''} left"
        parts.append(s + ".")
    return " ".join(parts)


def symbol_svg(sym):
    """Lewis symbol: one dot per side first, then pairs (top, right, bottom, left)."""
    v = VALENCE[sym]
    side = [0, 0, 0, 0]
    for i in range(v):
        side[i % 4] += 1
    color = _var(sym)
    angles = [-90, 0, 90, 180]
    dots = []
    for n, t in zip(side, angles):
        if n == 2:
            dots.append(_dot_pair(0, 0, t, color, sym))
        elif n == 1:
            dots.append(_dot_single(0, 0, t, color, sym))
    pairs = sum(1 for n in side if n == 2)
    singles = sum(1 for n in side if n == 1)
    alt = (f"Lewis symbol for {WORD[sym]}: {v} valence electron{'s' if v > 1 else ''}, "
           f"{pairs} pair{'s' if pairs != 1 else ''} and {singles} single electron{'s' if singles != 1 else ''}")
    return (f'<svg class="lewis sym" viewBox="-30 -26 60 52" width="66" height="58" role="img" '
            f'aria-labelledby="t-sym-{sym}" xmlns="http://www.w3.org/2000/svg"><title id="t-sym-{sym}">{alt}</title>'
            f'{"".join(dots)}<g class="atoms"><text x="0" y="0" dy="0.35em" style="fill:{color}">{sym}</text></g></svg>'), pairs, singles
