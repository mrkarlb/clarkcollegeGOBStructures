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
    if problems:
        print("STRUCTURE CHECK FAILED:\n" + "\n".join(problems))
        sys.exit(1)
    print(f"Structure check passed: {len(S)} structures ({len(named)} name-checked), "
          f"{len(PRACTICE)} practice problems.")


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


# ---------------------------------------------------------------- markdown
def expand_shortcodes(md):
    def repl(m):
        kind, *args = m.group(1).split()
        if kind == "lewis":
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
PAGE = PAGES[0]


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
