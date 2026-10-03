# Molecular Structure (CHEM&121)

A two-page web guide for Clark College CHEM&121 students preparing for health-profession careers:

- **Lewis Structures** (`index.html`): Lewis symbols, drawing, bonding patterns and formal charge, resonance
- **Shape, Polarity, and Intermolecular Forces** (`shape.html`): VSEPR, bond and molecular polarity,
  intermolecular forces, beyond the octet

Both pages open with the same roadmap (Lewis structure → shape → bond polarity → molecular polarity →
intermolecular forces → properties) and link to each other. CHEM&131 students use Bonding Patterns,
Resonance, and Intermolecular Forces as review.

Written by Dr. Karl Bailey, Clark College. Licensed CC BY-NC 4.0 (see `LICENSE`).

## How it works

| File | What it holds |
|---|---|
| `content/lewis/*.md`, `content/shape/*.md` | The text of each page, one Markdown file per section, in page order |
| `structures.py` | Every structure: id, formula, name, SMILES with formal charges, and layout; plus the practice problems |
| `lewis.py` | Works out lone pairs, checks electron counts, and draws each Lewis structure as an SVG that follows the page's light/dark theme |
| `shapes.py` | Works out VSEPR shape, partial charges, and polarity from the same structures; draws 3D shapes and the intermolecular-force diagrams |
| `build.py` | Checks every structure, builds both pages into `_site/`, and checks every link |
| `template.html` | Page layout, styles, search, and dark mode |
| `tools/review_sheet.py` | Builds `_site/structure_review.html`, every structure on one page with its electron counts, for chemistry review |

Every push to `main` runs `.github/workflows/pages.yml`, which builds and publishes the site.

**The build fails if any structure is wrong.** Each structure is entered as SMILES with formal
charges, and the build works out every lone pair from that. It then checks that:

- the atoms match the formula, and the formal charges add up to the charge;
- every H has 2 electrons, and every C, N, O, and F has exactly 8;
- S, P, Xe and other larger atoms have at least 8 electrons in no more than six groups;
- no atom has an unpaired electron;
- the electrons drawn equal the valence electrons counted from the formula;
- the name, read by OPSIN, gives the same formula and charge;
- RDKit didn't have to rewrite the structure to read it;
- any molecule the text calls polar or nonpolar (`polar=`) gets that answer from the course rules,
  unless `polar_note=` records a real exception (ozone).

Shape names, bond angles, ΔEN values, and partial charges on the page are all computed from the
structures, so they can't disagree with the drawings. The build also fails if any link to a section
or to the other page points nowhere.

Electron totals in the text come from `{{ve id}}`, computed from the formula, so they can't
disagree with the drawings.

## Editing

- **Change text:** edit the Markdown file in `content/`. Structures go in with shortcodes on their own line:
  - `[[lewis id]]` one structure; `[[lewisrow id1 id2]]` a row
  - `[[steps id]]` Steps 1–3 for a neutral molecule (add `share` when a lone pair has to be shared)
  - `[[resonance id1 id2 …]]` resonance forms joined by ↔
  - `[[count id]]` the valence-electron count table; `[[symbols C N O]]` Lewis symbols
  - `[[patterns]]` the bonding-pattern table; `[[practice id]]` a practice problem
  - `[[roadmap]]` the topic roadmap, with this page's steps highlighted
  - `[[shapetable id…]]`, `[[shape id]]`, `[[shaperow id…]]`, `[[shapeq id]]` VSEPR table, 3D shapes, shape practice
  - `[[en]]`, `[[bondtable C-H C=O …]]`, `[[polar id…]]` electronegativity, bond polarity, molecular polarity
  - `[[hydration]]`, `[[hbond]]`, `[[extshapes]]` ion–dipole and hydrogen-bond diagrams; five- and six-group shapes
- **Link to the other page:** `shape.html#imf` or `index.html#patterns`; section anchors come from the `{#id}` on each heading.
- **Add a structure:** add an `L(...)` line to `structures.py`. Resonance forms use `like=` to share
  the first form's layout, with atoms written in the same order.
- **Add a practice problem:** add a line to `PRACTICE`.

## Building locally (optional)

Requires Python 3.11+ and Java.

```
pip install -r requirements.txt
python build.py
python tools/review_sheet.py
```

Then open `_site/index.html` or `_site/shape.html`.
