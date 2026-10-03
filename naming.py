"""Ions, compounds, and acids for the Naming Compounds page.

Nothing on that page is typed by hand twice. Each ion is entered once with its
charge; every ionic formula is worked out from the charges, every name from the
naming rules, and every acid from its anion. The build then checks each one:

  - ionic compounds: charges balance, and OPSIN reads the generated name as the
    same formula (and, for each chart ion, as the same ion and charge)
  - covalent compounds: the prefixes in the name match the atom counts in the
    formula, and OPSIN agrees whenever it can read the name
  - acids: OPSIN reads the acid name as H⁺ plus the anion

Any failure stops the build.
"""
import math
import re
from fractions import Fraction

from lewis import formula_html

# ---------------------------------------------------------------- ions
# kind: "main" (charge from the group), "fixed" (transition metal with one common charge),
#       "variable" (needs a Roman numeral), "poly" (polyatomic)
# paren: polyatomic ions get parentheses when more than one is needed.
IONS = []


def ion(id, formula, name, kind, old="", aka="", where="", group=None, paren=None, family=""):
    body, _, chg = formula.partition("^")
    m = re.fullmatch(r"(\d*)([+-])", chg)
    charge = (int(m.group(1)) if m.group(1) else 1) * (1 if m.group(2) == "+" else -1)
    atoms = sum(int(n or 1) for _, n in re.findall(r"([A-Z][a-z]?)(\d*)", body))
    IONS.append(dict(id=id, formula=formula, body=body, charge=charge, name=name, kind=kind, old=old, aka=aka,
                     where=where, group=group, family=family,
                     paren=(kind == "poly" and atoms > 1) if paren is None else paren))


# main-group cations: the ion takes the element's name
ion("li", "Li^+", "lithium", "main", group=1, where="lithium carbonate, a mood-stabilizing medicine")
ion("na", "Na^+", "sodium", "main", group=1, where="the main positive ion outside your cells")
ion("k", "K^+", "potassium", "main", group=1, where="the main positive ion inside your cells")
ion("mg", "Mg^2+", "magnesium", "main", group=2)
ion("ca", "Ca^2+", "calcium", "main", group=2, where="bones, teeth, muscle contraction")
ion("ba", "Ba^2+", "barium", "main", group=2)
ion("al", "Al^3+", "aluminum", "main", group=13)
# main-group anions: element name with -ide
ion("n3", "N^3-", "nitride", "main", group=15)
ion("p3", "P^3-", "phosphide", "main", group=15)
ion("o2", "O^2-", "oxide", "main", group=16)
ion("s2", "S^2-", "sulfide", "main", group=16)
ion("f", "F^-", "fluoride", "main", group=17, where="toothpaste and drinking water")
ion("cl", "Cl^-", "chloride", "main", group=17, where="the main negative ion outside your cells")
ion("br", "Br^-", "bromide", "main", group=17)
ion("i", "I^-", "iodide", "main", group=17, where="iodized salt; your thyroid uses it")
# transition metals with only one common charge: no Roman numeral
ion("ag", "Ag^+", "silver", "fixed", where="silver-containing wound dressings")
ion("zn", "Zn^2+", "zinc", "fixed", where="zinc oxide sunscreen and diaper cream")
# transition (and post-transition) metals with more than one charge: Roman numeral required
ion("fe2", "Fe^2+", "iron(II)", "variable", old="ferrous", where="iron supplements (ferrous sulfate)")
ion("fe3", "Fe^3+", "iron(III)", "variable", old="ferric")
ion("cu1", "Cu^+", "copper(I)", "variable", old="cuprous")
ion("cu2", "Cu^2+", "copper(II)", "variable", old="cupric")
ion("sn2", "Sn^2+", "tin(II)", "variable", old="stannous", where="stannous fluoride toothpaste")
ion("sn4", "Sn^4+", "tin(IV)", "variable", old="stannic")
ion("pb2", "Pb^2+", "lead(II)", "variable", old="plumbous")
ion("pb4", "Pb^4+", "lead(IV)", "variable", old="plumbic")
ion("cr2", "Cr^2+", "chromium(II)", "variable", old="chromous")
ion("cr3", "Cr^3+", "chromium(III)", "variable", old="chromic")
ion("hg1", "Hg2^2+", "mercury(I)", "variable", old="mercurous", paren=False,
    aka="two Hg atoms bonded together, each +1")
ion("hg2", "Hg^2+", "mercury(II)", "variable", old="mercuric")
# polyatomic ions
ion("nh4", "NH4^+", "ammonium", "poly", family="cation", where="ammonium salts in smelling salts and fertilizers")
ion("h3o", "H3O^+", "hydronium", "poly", family="cation", where="what makes a solution acidic")
ion("oh", "OH^-", "hydroxide", "poly", family="other", where="antacids such as milk of magnesia")
ion("cn", "CN^-", "cyanide", "poly", family="other", where="a poison that blocks your cells from using oxygen")
ion("ac", "CH3COO^-", "acetate", "poly", family="other", aka="also written C₂H₃O₂⁻", where="vinegar")
ion("o22", "O2^2-", "peroxide", "poly", family="other", paren=False, where="hydrogen peroxide disinfectant")
ion("c2o4", "C2O4^2-", "oxalate", "poly", family="other", where="calcium oxalate kidney stones")
ion("no3", "NO3^-", "nitrate", "poly", family="N")
ion("no2", "NO2^-", "nitrite", "poly", family="N", where="a preservative in cured meats")
ion("co3", "CO3^2-", "carbonate", "poly", family="C", where="antacids and bones")
ion("hco3", "HCO3^-", "hydrogen carbonate", "poly", family="C", aka="bicarbonate",
    where="carries CO₂ in your blood and buffers its pH")
ion("so4", "SO4^2-", "sulfate", "poly", family="S")
ion("hso4", "HSO4^-", "hydrogen sulfate", "poly", family="S", aka="bisulfate")
ion("so3", "SO3^2-", "sulfite", "poly", family="S", where="a preservative in wine and dried fruit")
ion("hso3", "HSO3^-", "hydrogen sulfite", "poly", family="S", aka="bisulfite")
ion("po4", "PO4^3-", "phosphate", "poly", family="P", where="bones, teeth, DNA, and ATP")
ion("hpo4", "HPO4^2-", "hydrogen phosphate", "poly", family="P", aka="monohydrogen phosphate",
    where="buffers the pH inside your cells")
ion("h2po4", "H2PO4^-", "dihydrogen phosphate", "poly", family="P", where="buffers the pH inside your cells")
ion("clo", "ClO^-", "hypochlorite", "poly", family="Cl", where="household bleach")
ion("clo2", "ClO2^-", "chlorite", "poly", family="Cl")
ion("clo3", "ClO3^-", "chlorate", "poly", family="Cl")
ion("clo4", "ClO4^-", "perchlorate", "poly", family="Cl")

ION = {i["id"]: i for i in IONS}


def ion_html(i):
    return formula_html(i["formula"])


def _count(n):
    return "" if n == 1 else str(n)


def ionic(cat, an):
    """Formula, name, and counts for the compound of two ions; charges balance by construction."""
    c, a = ION[cat], ION[an]
    if c["charge"] <= 0 or a["charge"] >= 0:
        raise ValueError(f"{cat}+{an}: need a cation and an anion")
    lcm = c["charge"] * -a["charge"] // math.gcd(c["charge"], -a["charge"])
    nc, na = lcm // c["charge"], lcm // -a["charge"]
    part = lambda i, n: (f"({i['body']}){n}" if i["paren"] and n > 1 else i["body"] + _count(n))
    return dict(id=f"{cat}-{an}", cat=c, an=a, nc=nc, na=na, lcm=lcm,
                formula=part(c, nc) + part(a, na), name=f"{c['name']} {a['name']}")


# ---------------------------------------------------------------- formulas
def parse(formula):
    """'Ca3(PO4)2' -> {'Ca': 3, 'P': 2, 'O': 8}"""
    def walk(s, i=0):
        counts = {}
        while i < len(s):
            if s[i] == "(":
                inner, i = walk(s, i + 1)
                n = re.match(r"\d*", s[i:]).group()
                i += len(n)
                for k, v in inner.items():
                    counts[k] = counts.get(k, 0) + v * int(n or 1)
            elif s[i] == ")":
                return counts, i + 1
            else:
                m = re.match(r"([A-Z][a-z]?)(\d*)", s[i:])
                if not m:
                    raise ValueError(f"can't read formula {formula}")
                counts[m.group(1)] = counts.get(m.group(1), 0) + int(m.group(2) or 1)
                i += len(m.group())
        return counts, i
    return walk(formula)[0]


def html_formula(f):
    return formula_html(f)


# ---------------------------------------------------------------- covalent
PREFIX = ["mono", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca"]
ELEMENT = {"hydrogen": "H", "boron": "B", "carbon": "C", "nitrogen": "N", "oxygen": "O", "fluorine": "F",
           "silicon": "Si", "phosphorus": "P", "sulfur": "S", "chlorine": "Cl", "selenium": "Se",
           "bromine": "Br", "iodine": "I", "xenon": "Xe"}
IDE = {"hydride": "H", "carbide": "C", "nitride": "N", "oxide": "O", "fluoride": "F", "phosphide": "P",
       "sulfide": "S", "chloride": "Cl", "selenide": "Se", "bromide": "Br", "iodide": "I"}

COVALENT = []


def cov(id, formula, name, common="", where=""):
    COVALENT.append(dict(id=id, formula=formula, name=name, common=common, where=where))


cov("co", "CO", "carbon monoxide", where="a poison that blocks hemoglobin")
cov("co2", "CO2", "carbon dioxide", where="what you exhale")
cov("n2o", "N2O", "dinitrogen monoxide", common="nitrous oxide", where="laughing gas, a dental anesthetic")
cov("no", "NO", "nitrogen monoxide", common="nitric oxide", where="a signal that relaxes blood vessels")
cov("no2", "NO2", "nitrogen dioxide", where="smog")
cov("n2o4", "N2O4", "dinitrogen tetroxide")
cov("so2", "SO2", "sulfur dioxide", where="a preservative in dried fruit")
cov("so3", "SO3", "sulfur trioxide")
cov("p2o5", "P2O5", "diphosphorus pentoxide")
cov("pcl3", "PCl3", "phosphorus trichloride")
cov("pcl5", "PCl5", "phosphorus pentachloride")
cov("ccl4", "CCl4", "carbon tetrachloride")
cov("sf6", "SF6", "sulfur hexafluoride", where="a greenhouse gas")
cov("nf3", "NF3", "nitrogen trifluoride")
cov("cs2", "CS2", "carbon disulfide")
cov("of2", "OF2", "oxygen difluoride")
COV = {c["id"]: c for c in COVALENT}

COMMON = [("H2O", "water"), ("NH3", "ammonia"), ("CH4", "methane"), ("H2O2", "hydrogen peroxide")]


def read_prefixed(name):
    """'dinitrogen tetroxide' -> [('N', 2), ('O', 4)], following the textbook rules; raises if a rule is broken."""
    first, second = name.split()

    def split(word, ends):
        for k in sorted(range(len(PREFIX)), key=lambda k: -len(PREFIX[k])):
            p = PREFIX[k]
            if word.startswith(p) and word[len(p):] in ends:
                return k + 1, word[len(p):], p
            if p[-1] in "ao" and word.startswith(p[:-1]) and word[len(p) - 1:] in ends and word[len(p) - 1] in "ao":
                return k + 1, word[len(p) - 1:], p[:-1]
        return (1, word, "") if word in ends else (None, word, "")

    n1, root1, p1 = split(first, ELEMENT)
    n2, root2, p2 = split(second, IDE)
    if n1 is None or n2 is None:
        raise ValueError(f"can't read '{name}'")
    if p1 == "mono":
        raise ValueError(f"'{name}': the first element never takes mono-")
    if n2 == 1 and not p2:
        raise ValueError(f"'{name}': the second element always takes a prefix, even mono-")
    if p2 in PREFIX and p2[-1] in "ao" and root2[0] in "aeiou" and root2 == "oxide":
        raise ValueError(f"'{name}': drop the prefix's final vowel before oxide ({p2[:-1]}{root2})")
    return [(ELEMENT[root1], n1), (IDE[root2], n2)]


def prefixed_name(formula):
    """Write the prefix name from a two-element formula, the way students are taught to."""
    (e1, n1), (e2, n2) = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    n1, n2 = int(n1 or 1), int(n2 or 1)
    w1 = next(k for k, v in ELEMENT.items() if v == e1)
    w2 = next(k for k, v in IDE.items() if v == e2)
    p1 = "" if n1 == 1 else PREFIX[n1 - 1]
    p2 = PREFIX[n2 - 1]
    if p2[-1] in "ao" and w2 == "oxide":
        p2 = p2[:-1]
    return f"{p1}{w1} {p2}{w2}"


# ---------------------------------------------------------------- acids
ACID_ROOT = {"sulf": "sulfur", "phosph": "phosphor"}


def acid(an):
    a = ION[an]
    nm = a["name"]
    if nm.endswith("ide"):
        root = nm[:-3]
        name = f"hydro{ACID_ROOT.get(root, root)}ic acid"
        rule = "-ide becomes hydro- … -ic acid"
    elif nm.endswith("ate"):
        root = nm[:-3]
        name = f"{ACID_ROOT.get(root, root)}ic acid"
        rule = "-ate becomes -ic acid"
    elif nm.endswith("ite"):
        root = nm[:-3]
        name = f"{ACID_ROOT.get(root, root)}ous acid"
        rule = "-ite becomes -ous acid"
    else:
        raise ValueError(f"no acid rule for {nm}")
    n = -a["charge"]
    formula = "CH3COOH" if an == "ac" else f"H{_count(n)}{a['body']}"
    return dict(id=f"acid-{an}", an=a, name=name, formula=formula, rule=rule, nH=n)


# ---------------------------------------------------------------- verification
def _smiles_counts(smiles):
    from rdkit import Chem
    counts, charge = {}, 0
    for frag in smiles.split("."):
        m = Chem.MolFromSmiles(frag)
        m = Chem.AddHs(m)
        for at in m.GetAtoms():
            counts[at.GetSymbol()] = counts.get(at.GetSymbol(), 0) + 1
            charge += at.GetFormalCharge()
    return counts, charge


def _ratio(counts):
    g = 0
    for v in counts.values():
        g = math.gcd(g, v)
    return {k: v // g for k, v in counts.items()}


NO_OPSIN = {"p3": "checked by its group charge", "o22": "O₂ with a 2− charge; OPSIN doesn't read peroxide salts"}


def verify(compound_ids, acid_ids, cov_ids):
    """Returns a list of problems (empty when everything checks out)."""
    from py2opsin import py2opsin
    problems, queries = [], []

    # main-group charges follow the group: metals lose electrons (group 1, 2, 13), nonmetals gain (group − 18)
    for i in IONS:
        if i["kind"] == "main":
            want = i["group"] if i["group"] <= 2 else (3 if i["group"] == 13 else i["group"] - 18)
            if i["charge"] != want:
                problems.append(f"  ion {i['id']}: group {i['group']} gives {want:+d}, not {i['charge']:+d}")
    # every chart ion, through a simple compound OPSIN can read (OPSIN doesn't know these two names)
    for i in IONS:
        if i["id"] in NO_OPSIN:
            continue
        partner = ION["cl"] if i["charge"] > 0 else ION["na"]
        cmp = ionic(i["id"], "cl") if i["charge"] > 0 else ionic("na", i["id"])
        queries.append(("ion " + i["id"], cmp["name"], parse(cmp["formula"]), 0, True))
    for cid in compound_ids:
        cat, an = cid.split("-")
        cmp = ionic(cat, an)
        got = cmp["nc"] * cmp["cat"]["charge"] + cmp["na"] * cmp["an"]["charge"]
        if got:
            problems.append(f"  {cid}: charges add to {got}")
        queries.append((cid, cmp["name"], parse(cmp["formula"]), 0, True))
    for aid in acid_ids:
        a = acid(aid)
        want = parse(a["formula"])
        if parse(a["formula"]) != {**parse(a["an"]["body"]), "H": parse(a["an"]["body"]).get("H", 0) + a["nH"]}:
            problems.append(f"  acid {aid}: formula {a['formula']} isn't H⁺ plus {a['an']['formula']}")
        queries.append(("acid " + aid, a["name"], want, 0, False))
    for c in COVALENT:
        try:
            got = read_prefixed(c["name"])
        except ValueError as e:
            problems.append(f"  {c['id']}: {e}")
            continue
        want = [(e, int(n or 1)) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", c["formula"])]
        if got != want:
            problems.append(f"  {c['id']}: the name '{c['name']}' means {got}, not {c['formula']}")
        if prefixed_name(c["formula"]) != c["name"]:
            problems.append(f"  {c['id']}: the rules give '{prefixed_name(c['formula'])}', not '{c['name']}'")
        queries.append(("covalent " + c["id"], c["name"], parse(c["formula"]), 0, None))
    for cid in cov_ids:
        if cid not in COV:
            problems.append(f"  no covalent compound '{cid}'")

    results = py2opsin([q[1] for q in queries])
    read = 0
    for (label, name, want, charge, ratio), smi in zip(queries, results):
        if not smi:
            if ratio is None:      # covalent: the prefix check above is the authority
                continue
            problems.append(f"  {label}: OPSIN could not read '{name}'")
            continue
        got, q = _smiles_counts(smi)
        same = (_ratio(got) == _ratio(want)) if ratio else (got == want)
        if not same or q != charge:
            problems.append(f"  {label}: OPSIN reads '{name}' as {smi}, which is not {want}")
        read += 1
    return problems, read, len(queries)
