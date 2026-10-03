"""Build the CHEM&121 molecular structure site.

    python build.py            # check every structure, render, write _site/index.html

The build stops if any structure fails its electron bookkeeping or its name
check, so a wrong structure can never be published.
"""
import datetime
import glob
import html
import os
import re
import sys
import warnings
from zoneinfo import ZoneInfo

import markdown
from rdkit import Chem, RDLogger

import lewis
import naming
import shapes
from lewis import (audit, charge_text, describe, electron_count_rows, formula_html, lewis_svg,
                   make_mol, parse_formula, symbol_svg)
from structures import PRACTICE, S

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "_site")

BY_ID = {s["id"]: s for s in S}
PRACTICE_BY_ID = {p[0]: p for p in PRACTICE}
MOLS = {}


# ---------------------------------------------------------------- verification
def _formula_of(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    counts = {}
    for a in m.GetAtoms():
        counts[a.GetSymbol()] = counts.get(a.GetSymbol(), 0) + 1
    return counts, sum(a.GetFormalCharge() for a in m.GetAtoms())


def verify():
    problems = []
    if len(BY_ID) != len(S):
        problems.append("  duplicate structure ids")
    for s in S:  # lay out every structure (resonance forms borrow their partner's layout)
        try:
            like = MOLS[s["like"]] if s["like"] else None
            MOLS[s["id"]] = make_mol(s["smiles"], s["coords"], like, s["rotate"], s["flip"])
        except (ValueError, KeyError) as e:
            problems.append(f"  {s['id']}: {e}")
            continue
        if s["formula"] or s["role"] == "pattern":
            problems += ["  " + p for p in audit_any(s)]
    named = [s for s in S if s["opsin"]]
    if named:
        from py2opsin import py2opsin
        for s, o in zip(named, py2opsin([s["opsin"] for s in named])):
            if not o:
                problems.append(f"  {s['id']}: OPSIN could not read the name '{s['opsin']}'")
            elif _formula_of(o) != parse_formula(s["formula"]):
                problems.append(f"  {s['id']}: OPSIN reads '{s['opsin']}' as {o}, which is not {s['formula']}")
    for s in S:
        if s["polar"] is None or s["id"] not in MOLS:
            continue
        got, why = shapes.polarity(MOLS[s["id"]])
        want = "polar" if s["polar"] else "nonpolar"
        if got != want and not s["polar_note"]:
            problems.append(f"  {s['id']}: marked {want}, but the course rules give {got} ({why})")
        if got == want and s["polar_note"]:
            problems.append(f"  {s['id']}: polar_note set, but the rules already give {got}; remove the note")
    for pid, ids, _ in PRACTICE:
        for i in ids:
            if i not in BY_ID:
                problems.append(f"  practice {pid}: no structure '{i}'")
    nprob, nread, nall = verify_naming()
    problems += nprob
    if problems:
        print("STRUCTURE CHECK FAILED:\n" + "\n".join(problems))
        sys.exit(1)
    print(f"Structure check passed: {len(S)} structures ({len(named)} name-checked), "
          f"{len(PRACTICE)} practice problems.")
    print(f"Naming check passed: {len(naming.IONS)} ions; {nall} names checked, {nread} also read by OPSIN.")


def verify_naming():
    """Collect every compound and acid the naming page shows, then check them all."""
    ionic_ids, acid_ids, cov_ids = set(), set(), set()
    for path in glob.glob(os.path.join(ROOT, "content", "*", "*.md")):
        for kind, args in re.findall(r"^\[\[(ionic|acids|covalent|nameq|formulaq)\s+([^\]]+)\]\]", open(path).read(), re.M):
            for a in args.split():
                if kind == "acids" or a.startswith("acid-"):
                    acid_ids.add(a.removeprefix("acid-"))
                elif kind == "covalent" or a in naming.COV:
                    cov_ids.add(a)
                else:
                    ionic_ids.add(a)
    problems, read, total = naming.verify(sorted(ionic_ids), sorted(acid_ids), sorted(cov_ids))
    return problems, read, total


def audit_any(s):
    m = MOLS[s["id"]]
    if s["formula"]:
        return audit(m, s["formula"], s["id"])
    # pattern fragment: check only the atoms that are drawn
    probs = []
    for a in m.GetAtoms():
        if not a.GetAtomicNum():
            continue
        shell = 2 * lewis._bond_sum(a) + 2 * lewis._lone_pairs(a)
        if shell != 8:
            probs.append(f"{s['id']}: {a.GetSymbol()} has {shell} electrons")
    return probs


# ---------------------------------------------------------------- pieces
USED, SHOWN = set(), {}


def _uid(sid, suffix=""):
    SHOWN[sid + suffix] = SHOWN.get(sid + suffix, 0) + 1
    n = SHOWN[sid + suffix]
    return f"{sid}{suffix}" + ("" if n == 1 else f"-{n}")


def caption(s, show_note=True):
    f = f'<span class="fm">{formula_html(s["formula"])}</span> ' if s["formula"] else ""
    note = f'<span class="cm">{html.escape(s["note"])}</span>' if s["note"] and show_note else ""
    tag = '<span class="tag health">health</span>' if s["role"] == "health" else ""
    return f'<figcaption><span class="nm">{f}{html.escape(s["title"])}</span>{note}{tag}</figcaption>'


def drawing(sid, stage="final", partial=False):
    USED.add(sid)
    s = BY_ID[sid]
    name = f"{s['formula']} ({s['title']})" if s["formula"] else s["title"]
    name = name.replace("^", " ")
    notes = shapes.partials(MOLS[sid]) if partial else None
    uid = _uid(sid, ("" if stage == "final" else "-" + stage) + ("-partial" if partial else ""))
    return lewis_svg(uid, MOLS[sid], describe(MOLS[sid], name, stage, notes), stage,
                     brackets=s["role"] != "pattern", notes=notes)


def figure(sid, partial=False):
    return f'<figure class="fig"><div class="pic">{drawing(sid, partial=partial)}</div>{caption(BY_ID[sid])}</figure>'


# ---------------------------------------------------------------- page 2 pieces
def shape_fig(sid):
    USED.add(sid)
    s = BY_ID[sid]
    svg, g = shapes.shape_svg(_uid(sid, "-shape"), MOLS[sid], f"{s['formula'].replace('^', ' ')} ({s['title']})")
    lp = f"{g['lone_pairs']} lone pair{'s' if g['lone_pairs'] != 1 else ''}"
    return (f'<figure class="fig"><div class="pic">{svg}</div><figcaption><span class="nm">'
            f'<span class="fm">{formula_html(s["formula"])}</span> {html.escape(s["title"])}</span>'
            f'<span class="cm"><strong>{g["molecular"]}</strong> · {g["groups"]} groups, {lp}</span></figcaption></figure>')


def polar_fig(sid):
    USED.add(sid)
    s = BY_ID[sid]
    m = MOLS[sid]
    notes = shapes.partials(m)
    name = f"{s['formula']} ({s['title']})".replace("^", " ")
    svg = lewis_svg(_uid(sid, "-polar"), m, describe(m, name, "final", notes), notes=notes)
    got, why = shapes.polarity(m)
    verdict = got if not s["polar_note"] else ("polar" if s["polar"] else "nonpolar")
    reason = s["polar_note"] or why
    return (f'<figure class="fig"><div class="pic">{svg}</div><figcaption><span class="nm">'
            f'<span class="fm">{formula_html(s["formula"])}</span> {html.escape(s["title"])}</span>'
            f'<span class="verdict {verdict}">{verdict}</span><span class="cm">{html.escape(reason)}</span>'
            f'</figcaption></figure>')


def shape_table(ids):
    rows = []
    for sid in ids:
        USED.add(sid)
        s = BY_ID[sid]
        svg, g = shapes.shape_svg(_uid(sid, "-shape"), MOLS[sid], f"{s['formula'].replace('^', ' ')} ({s['title']})")
        rows.append(f'<tr><td>{g["groups"]}</td><td>{g["lone_pairs"]}</td><td>{g["electron"]}</td>'
                    f'<th scope="row">{g["molecular"]}</th><td>{g["angle"]}</td>'
                    f'<td class="shape-ex"><div class="pic">{svg}</div><span class="fm">{formula_html(s["formula"])}</span> '
                    f'{html.escape(s["title"])}</td></tr>')
    return ('<div class="table-wrap"><table class="shapes"><caption>Shapes for up to four electron groups. '
            "Count every atom bonded to the central atom (single, double, or triple bond: each counts once) "
            'and every lone pair on it.</caption><thead><tr><th scope="col">Electron groups</th>'
            '<th scope="col">Lone pairs</th><th scope="col">Electron geometry</th><th scope="col">Molecular geometry</th>'
            '<th scope="col">Bond angle</th><th scope="col">Example</th></tr></thead><tbody>'
            + "".join(rows) + "</tbody></table></div>")


def ext_shape_table():
    ex = {(5, 0): "PCl₅", (5, 1): "SF₄", (5, 2): "ClF₃", (5, 3): "XeF₂",
          (6, 0): "SF₆", (6, 1): "XeOF₄, BrF₅", (6, 2): "XeF₄"}
    rows = "".join(f'<tr><td>{k[0]}</td><td>{k[1]}</td><td>{shapes.GEOMETRY[k][0]}</td>'
                   f'<th scope="row">{shapes.GEOMETRY[k][1]}</th><td>{shapes.GEOMETRY[k][2]}</td><td>{v}</td></tr>'
                   for k, v in ex.items())
    return ('<div class="table-wrap"><table class="shapes ext"><caption>Beyond this course: five and six electron '
            'groups, found only around larger central atoms such as P, S, Cl, and Xe.</caption><thead><tr>'
            '<th scope="col">Electron groups</th><th scope="col">Lone pairs</th><th scope="col">Electron geometry</th>'
            '<th scope="col">Molecular geometry</th><th scope="col">Bond angles</th><th scope="col">Examples</th>'
            f"</tr></thead><tbody>{rows}</tbody></table></div>")


def en_table():
    els = ["H", "C", "N", "O", "F", "P", "S", "Cl"]
    cells = "".join(f'<td><span class="el" style="color:var({lewis.COLOR[e]})">{e}</span><br>{shapes.EN[e]}</td>' for e in els)
    return ('<div class="table-wrap"><table class="en"><caption>Electronegativity values used in this course '
            '(Pauling scale)</caption><tbody><tr>' + cells + "</tr></tbody></table></div>")


def bond_table(bonds):
    rows = []
    for spec in bonds:
        a, sep, c = re.match(r"([A-Z][a-z]?)([-=#])([A-Z][a-z]?)", spec).groups()
        d = abs(shapes.EN[a] - shapes.EN[c])
        kind = shapes.bond_type(a, c)
        label = {"nonpolar": "nonpolar covalent", "polar": "polar covalent", "ionic": "ionic"}[kind]
        if kind == "polar":
            neg = a if shapes.EN[a] > shapes.EN[c] else c
            pos = c if neg == a else a
            label += f" (δ− on {neg}, δ+ on {pos})"
        sym = {"-": "–", "=": "=", "#": "≡"}[sep]
        rows.append(f"<tr><th scope=\"row\">{a}{sym}{c}</th><td>{shapes.EN[a]} and {shapes.EN[c]}</td>"
                    f"<td>{d:.1f}</td><td>{label}</td></tr>")
    return ('<div class="table-wrap"><table class="bonds"><caption>Bond polarity from the difference in '
            'electronegativity (ΔEN)</caption><thead><tr><th scope="col">Bond</th><th scope="col">EN values</th>'
            '<th scope="col">ΔEN</th><th scope="col">Bond type</th></tr></thead><tbody>'
            + "".join(rows) + "</tbody></table></div>")


def shape_question(sid):
    USED.add(sid)
    s = BY_ID[sid]
    svg, g = shapes.shape_svg(_uid(sid, "-shapeq"), MOLS[sid], f"{s['formula'].replace('^', ' ')} ({s['title']})")
    lp = f"{g['lone_pairs']} lone pair{'s' if g['lone_pairs'] != 1 else ''}"
    same = "the same as the electron geometry" if g["electron"] == g["molecular"] else f"electron geometry {g['electron']}"
    return (f'<div class="practice"><p class="q"><strong>What shape is</strong> <span class="fm">{formula_html(s["formula"])}</span> '
            f'({html.escape(s["title"])})? Give the electron geometry and the molecular geometry.</p>'
            f'<details class="answer"><summary>Show answer</summary><div class="figrow">'
            f'<figure class="fig"><div class="pic">{drawing(sid)}</div><figcaption><span class="cm">Lewis structure</span></figcaption></figure>'
            f'<figure class="fig"><div class="pic">{svg}</div><figcaption><span class="cm">Shape</span></figcaption></figure></div>'
            f'<p><span class="ans">{g["molecular"]}</span> ({same}). The central {g["symbol"]} has '
            f'{g["groups"]} electron groups: {g["groups"] - g["lone_pairs"]} bonded atom{"s" if g["groups"] - g["lone_pairs"] != 1 else ""} '
            f'and {lp}. Bond angle: {g["angle"]}.</p></details></div>')


ROADMAP = [("index.html", "drawing", "Lewis structure"), ("shape.html", "shape", "Shape"),
           ("shape.html", "bond-polarity", "Bond polarity"), ("shape.html", "polarity", "Molecular polarity"),
           ("shape.html", "imf", "Intermolecular forces"), ("shape.html", "imf-properties", "Properties")]


def roadmap():
    items = []
    for k, (pg, anchor, label) in enumerate(ROADMAP):
        here = pg == PAGE["file"]
        href = f"#{anchor}" if here else f"{pg}#{anchor}"
        items.append(f'<li class="{"here" if here else ""}"><a href="{href}">{label}</a></li>')
    return ('<nav class="roadmap" aria-label="How the topics connect"><ol>' + "".join(items) + "</ol>"
            '<p class="cm">Each step uses the one before it. Highlighted steps are on this page.</p></nav>')


def steps(sid, share=False):
    s = BY_ID[sid]
    if parse_formula(s["formula"])[1]:
        raise ValueError(f"[[steps {sid}]]: step drawings are for neutral molecules")
    labels = [("skeleton", "Step 1", "Connect every atom with a single bond"),
              ("dots", "Step 2", "Add each atom's remaining electrons"),
              ("final", "Step 3", "Pair up single electrons; share a lone pair if an atom is still short"
               if share else "Pair up single electrons into bonds")]
    cells = []
    for stage, step, what in labels:
        cells.append(f'<figure class="fig step"><div class="pic">{drawing(sid, stage)}</div>'
                     f'<figcaption><span class="nm">{step}</span><span class="cm">{what}</span></figcaption></figure>')
    return f'<div class="figrow steps" role="group" aria-label="Steps for {html.escape(s["title"])}">{"".join(cells)}</div>'


def resonance(ids):
    s = BY_ID[ids[0]]
    parts = []
    for k, sid in enumerate(ids):
        if k:
            parts.append('<span class="res-arrow" aria-label="resonance with">↔</span>')
        parts.append(f'<div class="pic">{drawing(sid)}</div>')
    tag = '<span class="tag health">health</span>' if s["role"] == "health" else ""
    note = s["note"] if s["note"] and not s["note"].startswith("form") else ""
    note = f'<span class="cm">{html.escape(note)}</span>' if note else ""
    return (f'<figure class="fig res"><div class="res-row">{"".join(parts)}</div>'
            f'<figcaption><span class="nm"><span class="fm">{formula_html(s["formula"])}</span> '
            f'{html.escape(s["title"])}: {len(ids)} resonance forms</span>{note}{tag}</figcaption></figure>')


def count_table(sid):
    s = BY_ID[sid]
    rows, charge, total = electron_count_rows(s["formula"])
    body = "".join(f'<tr><td>{el}</td><td>{n}</td><td>{v}</td><td>{n * v}</td></tr>' for el, n, v, _ in rows)
    if charge:
        word = "add" if charge < 0 else "subtract"
        body += (f'<tr><td colspan="3">Charge {charge_text(charge)}: {word} {abs(charge)} '
                 f'electron{"s" if abs(charge) > 1 else ""}</td><td>{-charge:+d}</td></tr>')
    body += f'<tr class="total"><th scope="row" colspan="3">Valence electrons to place</th><td>{total}</td></tr>'
    return (f'<div class="table-wrap"><table class="count"><caption>Counting valence electrons: '
            f'{formula_html(s["formula"])}</caption><thead><tr><th scope="col">Element</th>'
            f'<th scope="col">Atoms</th><th scope="col">Valence e⁻ each</th><th scope="col">Total</th>'
            f'</tr></thead><tbody>{body}</tbody></table></div>')


def symbols(syms):
    cells = []
    for sym in syms:
        svg, pairs, singles = symbol_svg(sym)
        v = lewis.VALENCE[sym]
        detail = f"{singles} single{'s' if singles != 1 else ''}"
        if pairs:
            detail += f", {pairs} pair{'s' if pairs != 1 else ''}"
        cells.append(f'<figure class="fig sym"><div class="pic">{svg}</div><figcaption><span class="nm">'
                     f'{v} valence e⁻</span><span class="cm">{detail}</span></figcaption></figure>')
    return f'<div class="figrow">{"".join(cells)}</div>'


PATTERNS = [("C", "pat_c", None, "pat_cm"), ("N", "pat_n", "pat_np", "pat_nm"), ("O", "pat_o", "pat_op", "pat_om")]


def pattern_table():
    def cell(sid):
        if not sid:
            return '<td class="pat-none">not in this course</td>'
        m = MOLS[sid]
        a = next(x for x in m.GetAtoms() if x.GetAtomicNum())
        nb, lp = a.GetDegree(), lewis._lone_pairs(a)
        what = f"{nb} bond{'s' if nb != 1 else ''}, {lp} lone pair{'s' if lp != 1 else ''}"
        return f'<td><div class="pic">{drawing(sid)}</div><span class="cm">{what}</span></td>'
    rows = "".join(f'<tr><th scope="row">{el}</th>{cell(n)}{cell(p)}{cell(mn)}</tr>' for el, n, p, mn in PATTERNS)
    return ('<div class="table-wrap"><table class="patterns"><caption>Bonding patterns for C, N, and O. '
            "Count bonds, not electrons: a double bond counts as two bonds, a triple bond as three. "
            "Every pattern shown has 8 electrons.</caption><thead><tr><th scope=\"col\">Atom</th>"
            '<th scope="col">Neutral</th><th scope="col">+1 charge</th><th scope="col">−1 charge</th></tr></thead>'
            f"<tbody>{rows}</tbody></table></div>")


def practice(pid):
    _, ids, why = PRACTICE_BY_ID[pid]
    s = BY_ID[ids[0]]
    _, _, total = electron_count_rows(s["formula"])
    what = {1: "the Lewis structure of", 2: "both resonance forms of"}.get(len(ids), "all the resonance forms of")
    prompt = (f'<p class="q"><strong>Draw {what}</strong> <span class="fm">{formula_html(s["formula"])}</span> '
              f'({html.escape(s["title"])}).</p>')
    pics = resonance(ids) if len(ids) > 1 else f'<div class="pic">{drawing(ids[0])}</div>'
    answer = f'{pics}<p><strong>{total} valence electrons.</strong> {html.escape(why)}</p>'
    return (f'<div class="practice">{prompt}<details class="answer"><summary>Show answer</summary>'
            f"{answer}</details></div>")



# ---------------------------------------------------------------- naming page pieces
NM = naming


def _tbl(cls, caption, heads, rows):
    th = "".join(f'<th scope="col">{h}</th>' for h in heads)
    return (f'<div class="table-wrap"><table class="{cls}"><caption>{caption}</caption><thead><tr>{th}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def _ion(i):
    return f'<span class="fm">{NM.ion_html(i)}</span>'


def _q(q):
    """Charge with its number always shown: 1+, 2−."""
    return f"{abs(q)}{'+' if q > 0 else '−'}"


def _cap(t):
    return t[:1].upper() + t[1:]


NUMWORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}


def ion_chart(which):
    if which == "main":
        rows = []
        for g in (1, 2, 13, 15, 16, 17):
            ions = [i for i in NM.IONS if i["kind"] == "main" and i["group"] == g]
            q = ions[0]["charge"]
            kind = "metal: loses" if q > 0 else "nonmetal: gains"
            rows.append(f'<tr><th scope="row">{g}</th><td>{_q(q)}</td><td>{kind} {abs(q)} '
                        f'electron{"s" if abs(q) > 1 else ""}</td><td>'
                        + ", ".join(f'{_ion(i)} {i["name"]}' for i in ions) + "</td></tr>")
        return _tbl("ions", "Main-group ions: the charge comes from the group number", ["Group", "Charge", "Why", "Ions"], rows)
    if which == "variable":
        rows = [f'<tr><td>{_ion(i)}</td><th scope="row">{i["name"]}</th><td>{i["old"]}</td>'
                f'<td>{html.escape(i["where"] or i["aka"])}</td></tr>'
                for i in NM.IONS if i["kind"] == "variable"]
        return _tbl("ions", "Metals with more than one common charge: the Roman numeral gives the charge",
                    ["Ion", "Name", "Older name", "Notes"], rows)
    if which == "fixed":
        rows = [f'<tr><td>{_ion(i)}</td><th scope="row">{i["name"]}</th><td>{html.escape(i["where"])}</td></tr>'
                for i in NM.IONS if i["kind"] == "fixed"]
        return _tbl("ions", "Transition metals with only one common charge: no Roman numeral", ["Ion", "Name", "Where you'll see it"], rows)
    raise ValueError(f"[[ionchart {which}]]")


POLY_GROUPS = [("cation", "Positive polyatomic ions"), ("N", "Nitrogen"), ("C", "Carbon"), ("S", "Sulfur"),
               ("P", "Phosphorus"), ("Cl", "Chlorine"), ("other", "Others to know")]


def poly_chart():
    rows = []
    for fam, label in POLY_GROUPS:
        rows.append(f'<tr class="grp"><th scope="rowgroup" colspan="4">{label}</th></tr>')
        for i in NM.IONS:
            if i["kind"] == "poly" and i["family"] == fam:
                rows.append(f'<tr><td>{_ion(i)}</td><th scope="row">{i["name"]}</th>'
                            f'<td>{html.escape(i["aka"])}</td><td>{html.escape(i["where"])}</td></tr>')
    return _tbl("ions poly", "Common polyatomic ions. Use this chart; you don't need to memorize it.",
                ["Ion", "Name", "Also written or called", "Where you'll see it"], rows)


def oxy_series():
    ids = ["clo", "clo2", "clo3", "clo4"]
    pattern = {"clo": "hypo- … -ite", "clo2": "-ite", "clo3": "-ate", "clo4": "per- … -ate"}
    rows = [f'<tr><td>{_ion(NM.ION[i])}</td><td>{NM.parse(NM.ION[i]["body"])["O"]}</td>'
            f'<th scope="row">{NM.ION[i]["name"]}</th><td>{pattern[i]}</td></tr>' for i in ids]
    return _tbl("ions", "The chlorine series. Bromine and iodine follow the same pattern.",
                ["Ion", "O atoms", "Name", "Pattern"], rows)


def ionic_table(ids):
    rows = []
    for cid in ids:
        c = NM.ionic(*cid.split("-"))
        bal = (f'{c["nc"]} × ({_q(c["cat"]["charge"])}) = +{c["lcm"]}, '
               f'{c["na"]} × ({_q(c["an"]["charge"])}) = −{c["lcm"]}')
        rows.append(f'<tr><td>{_ion(c["cat"])} and {_ion(c["an"])}</td><td>{bal}</td>'
                    f'<td><span class="fm">{NM.html_formula(c["formula"])}</span></td><th scope="row">{c["name"]}</th></tr>')
    return _tbl("ionic", "Charges balance to zero", ["Ions", "Balance the charges", "Formula", "Name"], rows)


def prefix_table():
    cells = "".join(f'<tr><td>{n + 1}</td><th scope="row">{p}-</th></tr>' for n, p in enumerate(NM.PREFIX))
    return _tbl("prefixes", "Prefixes for covalent compounds", ["Number", "Prefix"], [cells])


def covalent_table(ids):
    rows = []
    for cid in ids:
        c = NM.COV[cid]
        notes = "; ".join(x for x in (f"also called {c['common']}" if c["common"] else "", c["where"]) if x)
        rows.append(f'<tr><td><span class="fm">{NM.html_formula(c["formula"])}</span></td>'
                    f'<th scope="row">{c["name"]}</th><td>{html.escape(notes)}</td></tr>')
    return _tbl("covalent", "Covalent compounds: prefixes give the number of each atom", ["Formula", "Name", "Notes"], rows)


def common_table():
    rows = [f'<tr><td><span class="fm">{NM.html_formula(f)}</span></td><th scope="row">{n}</th></tr>' for f, n in NM.COMMON]
    return _tbl("covalent", "Common names you'll use instead of prefix names", ["Formula", "Name"], rows)


def acid_table(ids):
    rows = []
    for aid in ids:
        a = NM.acid(aid)
        rows.append(f'<tr><td>{_ion(a["an"])} {a["an"]["name"]}</td><td>{a["rule"]}</td>'
                    f'<td><span class="fm">{NM.html_formula(a["formula"])}</span></td><th scope="row">{a["name"]}</th></tr>')
    return _tbl("acids", "Acids are named from their anions", ["Anion", "Rule", "Acid", "Name"], rows)


def _element_of(i):
    return re.sub(r"\(.*\)", "", i["name"])


def _explain_ionic(c, direction):
    cat, an = c["cat"], c["an"]
    if direction == "formula":
        why = f'{_cap(cat["name"])} is {_ion(cat)} and {an["name"]} is {_ion(an)}. '
        if c["nc"] == 1 and c["na"] == 1:
            why += "The charges are equal and opposite, so one of each balances. "
        else:
            why += (f'The charges balance at {c["lcm"]}: {c["nc"]} × ({_q(cat["charge"])}) = +{c["lcm"]} and '
                    f'{c["na"]} × ({_q(an["charge"])}) = −{c["lcm"]}. ')
        for i, n in ((cat, c["nc"]), (an, c["na"])):
            if i["paren"] and n > 1:
                why += f'{_cap(i["name"])} is a polyatomic ion, so it goes in parentheses with the {n} outside. '
        if cat["id"] == "hg1":
            why += "Mercury(I) is always Hg₂²⁺, two Hg atoms together, so its formula keeps the 2. "
        return why.strip()
    # formula -> name
    if cat["kind"] == "variable":
        el = _element_of(cat)
        why = (f'The anion is {an["name"]}, {_ion(an)}. ' +
               (f'{_cap(NUMWORD[c["na"]])} of them make −{c["lcm"]}, ' if c["na"] > 1 else f'It carries −{c["lcm"]}, ') +
               (f'so the {NUMWORD[c["nc"]]} {el} atoms make +{c["lcm"]}: each is {_q(cat["charge"])}. '
                if c["nc"] > 1 and cat["id"] != "hg1" else f'so {el} must be +{c["lcm"]}. ' if cat["id"] != "hg1"
                else "balanced by Hg₂²⁺, mercury(I). ") +
               f'{_cap(el)} has more than one common charge, so the name needs a Roman numeral: {c["name"]}.')
        return why
    if cat["kind"] == "poly":
        why = (f'{_cap(cat["name"])} is a polyatomic ion, {_ion(cat)}, and polyatomic ions never take a Roman numeral. '
               f'Name the cation, then the anion: {c["name"]}.')
    else:
        why = (f'{_cap(_element_of(cat))} has only one common charge, {_ion(cat)}, so there\'s '
               f'no Roman numeral. Name the cation, then the anion: {c["name"]}.')
    if c["nc"] > 1 or c["na"] > 1:
        why += " Ionic names never use prefixes: the charges already tell you how many of each."
    return why


def _explain_cov(c):
    (e1, n1), (e2, n2) = NM.read_prefixed(c["name"])
    first = "no prefix, because the first element never takes mono-" if n1 == 1 else f"{NM.PREFIX[n1 - 1]}-"
    why = (f'Two nonmetals, so it\'s covalent and the name uses prefixes. {n1} {e1}: {first}. '
           f'{n2} {e2}: {NM.PREFIX[n2 - 1]}-, and the second element ends in -ide. ')
    p = NM.PREFIX[n2 - 1]
    if p[-1] in "ao" and c["name"].split()[1].startswith(p[:-1] + "o"):
        why += f'Drop the final {p[-1]} before "oxide": {p[:-1]}oxide. '
    return why + (f'You\'ll also hear it called {c["common"]}.' if c["common"] else "")


def _explain_acid(a):
    nm = a["an"]["name"]
    extra = ""
    if nm.startswith(("sulf", "phosph")):
        extra = f' {_cap(nm[:-3])} adds a syllable in the acid: {a["name"].split()[0]}.'
    elif nm.startswith(("hypo", "per")):
        extra = f' The {"hypo-" if nm.startswith("hypo") else "per-"} prefix carries over to the acid.'
    return (f'The anion is {nm}, {_ion(a["an"])}. Rule: {a["rule"]}.{extra} '
            f'It takes {a["nH"]} H⁺ to balance the {_q(a["an"]["charge"])} charge.')


def naming_question(direction, ref):
    if ref.startswith("acid-"):
        a = NM.acid(ref[5:])
        formula, name, why = a["formula"], a["name"], _explain_acid(a)
    elif ref in NM.COV:
        c = NM.COV[ref]
        formula, name, why = c["formula"], c["name"], _explain_cov(c)
    else:
        c = NM.ionic(*ref.split("-"))
        formula, name, why = c["formula"], c["name"], _explain_ionic(c, direction)
    fm = f'<span class="fm">{NM.html_formula(formula)}</span>'
    if direction == "name":
        q, ans = f"<strong>Name</strong> {fm}.", f'<span class="ans">{name}</span>'
    else:
        q, ans = f"<strong>Write the formula for</strong> {name}.", f'<span class="ans">{NM.html_formula(formula)}</span>'
    return (f'<div class="practice"><p class="q">{q}</p><details class="answer"><summary>Show answer</summary>'
            f'<p>{ans}</p><p>{why}</p></details></div>')


NAMING_REFS = dict(ionic=set(), acid=set(), cov=set())

# ---------------------------------------------------------------- markdown
def expand_shortcodes(md):
    def repl(m):
        kind, *args = m.group(1).split()
        if kind == "ionchart":
            out = ion_chart(args[0])
        elif kind == "polychart":
            out = poly_chart()
        elif kind == "oxyseries":
            out = oxy_series()
        elif kind == "ionic":
            out = ionic_table(args)
        elif kind == "prefixes":
            out = prefix_table()
        elif kind == "covalent":
            out = covalent_table(args)
        elif kind == "commonnames":
            out = common_table()
        elif kind == "acids":
            out = acid_table(args)
        elif kind in ("nameq", "formulaq"):
            out = naming_question("name" if kind == "nameq" else "formula", args[0])
        elif kind == "lewis":
            out = figure(args[0])
        elif kind == "lewisrow":
            partial = "partial" in args
            out = '<div class="figrow">' + "".join(figure(a, partial) for a in args if a != "partial") + "</div>"
        elif kind == "steps":
            out = steps(args[0], share="share" in args[1:])
        elif kind == "resonance":
            out = resonance(args)
        elif kind == "count":
            out = count_table(args[0])
        elif kind == "symbols":
            out = symbols(args)
        elif kind == "patterns":
            out = pattern_table()
        elif kind == "practice":
            out = practice(args[0])
        elif kind == "shape":
            out = shape_fig(args[0])
        elif kind == "shaperow":
            out = '<div class="figrow">' + "".join(shape_fig(a) for a in args) + "</div>"
        elif kind == "shapetable":
            out = shape_table(args)
        elif kind == "extshapes":
            out = ext_shape_table()
        elif kind == "polar":
            out = '<div class="figrow">' + "".join(polar_fig(a) for a in args) + "</div>"
        elif kind == "en":
            out = en_table()
        elif kind == "roadmap":
            out = roadmap()
        elif kind == "hydration":
            out = (f'<figure class="fig wide"><div class="pic">{shapes.hydration_svg()}</div><figcaption>'
                   '<span class="nm">Ion–dipole attractions</span><span class="cm">Water surrounds each ion, '
                   'turning its oppositely charged end toward it. Lone pairs are left off for clarity.</span></figcaption></figure>')
        elif kind == "hbond":
            out = (f'<figure class="fig"><div class="pic">{shapes.hbond_svg()}</div><figcaption>'
                   '<span class="nm">Hydrogen bonds in water</span><span class="cm">Dotted lines are hydrogen bonds: '
                   'attractions between molecules, not covalent bonds.</span></figcaption></figure>')
        elif kind == "bondtable":
            out = bond_table(args)
        elif kind == "shapeq":
            out = shape_question(args[0])
        else:
            raise ValueError(f"unknown shortcode [[{m.group(1)}]]")
        return f"\n\n{out}\n\n"

    md = re.sub(r"^\[\[(.+?)\]\]\s*$", repl, md, flags=re.M)
    # inline: {{ve id}} = valence electron total, computed from the formula
    return re.sub(r"\{\{ve (\w+)\}\}", lambda m: str(electron_count_rows(BY_ID[m.group(1)]["formula"])[2]), md)


PAGES = [
    dict(file="naming.html", folder="naming", head="Naming Compounds — CHEM&amp;121 — Clark College",
         title="Naming Compounds", nav="Naming Compounds",
         sub="Ions, ionic and covalent compounds, and acids: how to go from a name to a formula and back. "
             "Written for students preparing for health-profession careers.",
         desc="Naming ions, ionic compounds, covalent compounds, and acids for students preparing for "
              "health-profession careers. CHEM&amp;121, Clark College."),
    dict(file="index.html", folder="lewis", head="Lewis Structures — CHEM&amp;121 — Clark College",
         title="Lewis Structures", nav="Lewis Structures",
         sub="Lewis symbols, drawing Lewis structures, formal charge, and resonance. Written for students "
             "preparing for health-profession careers.",
         desc="Lewis structures, formal charge, and resonance for students preparing for health-profession careers. "
              "CHEM&amp;121, Clark College."),
    dict(file="shape.html", folder="shape", head="Shape, Polarity, and Intermolecular Forces — CHEM&amp;121 — Clark College",
         title="Shape, Polarity, and Intermolecular Forces", nav="Shape, Polarity &amp; Forces",
         sub="From a Lewis structure to a molecule's shape, its polarity, and the forces between molecules. "
             "Written for students preparing for health-profession careers.",
         desc="Molecular shape (VSEPR), bond and molecular polarity, and intermolecular forces for students preparing "
              "for health-profession careers. CHEM&amp;121, Clark College."),
]
PAGE = PAGES[1]


def build_sections(folder):
    sections = []
    for path in sorted(glob.glob(os.path.join(ROOT, "content", folder, "*.md"))):
        md = open(path, encoding="utf-8").read()
        title, sid = re.match(r"#\s+(.+?)\s+\{#([\w-]+)\}", md).groups()
        body = expand_shortcodes(md)
        out = markdown.markdown(body, extensions=["tables", "attr_list", "md_in_html"])
        out = re.sub(r"<table>(.*?)</table>", r'<div class="table-wrap"><table>\1</table></div>', out, flags=re.S)
        sections.append((sid, title, out))
    return sections


def check_links(pages_html):
    """Every #anchor and page.html#anchor on the site must exist."""
    ids = {f: set(re.findall(r'\sid="([^"]+)"', h)) for f, h in pages_html.items()}
    bad = []
    for f, h in pages_html.items():
        for href in re.findall(r'href="([^"]+)"', h):
            if href.startswith("#"):
                if href[1:] not in ids[f]:
                    bad.append(f"  {f}: link to missing {href}")
            elif re.match(r"^[\w-]+\.html(#.*)?$", href):
                pg, _, anchor = href.partition("#")
                if pg not in ids:
                    bad.append(f"  {f}: link to missing page {pg}")
                elif anchor and anchor not in ids[pg]:
                    bad.append(f"  {f}: link to missing {href}")
    if bad:
        print("LINK CHECK FAILED:\n" + "\n".join(bad))
        sys.exit(1)
    print("Link check passed.")


def main():
    global PAGE
    verify()
    now = datetime.datetime.now(ZoneInfo("America/Los_Angeles"))
    template = open(os.path.join(ROOT, "template.html"), encoding="utf-8").read()
    os.makedirs(OUT, exist_ok=True)
    built = {}
    for PAGE in PAGES:
        sections = build_sections(PAGE["folder"])
        toc = "".join(f'<a href="#{sid}" data-target="{sid}">{html.escape(t)}</a>' for sid, t, _ in sections)
        main_html = "".join(f'<section class="block" id="{sid}-sec">{h}</section>' for sid, _, h in sections)
        nav = "".join(f'<a href="{p["file"]}"{" aria-current=page" if p is PAGE else ""}>{p["nav"]}</a>' for p in PAGES)
        page = template
        for key, val in [("{{HEAD}}", PAGE["head"]), ("{{DESC}}", PAGE["desc"]), ("{{TITLE}}", PAGE["title"]),
                         ("{{SUB}}", PAGE["sub"]), ("{{PAGENAV}}", nav), ("{{TOC}}", toc), ("{{MAIN}}", main_html),
                         ("{{UPDATED}}", now.strftime("%B %-d, %Y"))]:
            page = page.replace(key, val)
        built[PAGE["file"]] = page
        print(f"Built {PAGE['file']} ({len(page) // 1024} KB, {len(sections)} sections).")
    check_links(built)
    for f, page in built.items():
        with open(os.path.join(OUT, f), "w", encoding="utf-8") as fh:
            fh.write(page)
    open(os.path.join(OUT, ".nojekyll"), "w").close()
    unused = sorted(set(BY_ID) - USED)
    if unused:
        print("Note: structures not shown on either page:", ", ".join(unused))


if __name__ == "__main__":
    main()
