# Structure list for the CHEM&121 molecular structure site.
#
# Each structure is written as SMILES with formal charges. The build works out
# every lone pair from that, then checks each structure:
#   - the atoms match the formula, and formal charges add up to the charge
#   - every H has 2 electrons; every C, N, O, F has exactly 8
#   - S, P, Xe and other larger atoms have at least 8, in at most six groups
#   - no atom is left with an unpaired electron
#   - the electrons drawn equal the valence electrons counted from the formula
#   - OPSIN reads the name (opsin=...) as the same formula and charge
# Any failure stops the build, so a wrong structure can't be published.
#
# Fields
#   id       short id used in content shortcodes
#   formula  as students write it; charge after ^ (NH4^+, CO3^2-)
#   title    the name shown under the drawing
#   smiles   structure, with formal charges
#   note     short caption line (optional)
#   role     "worked", "health", "ref", "pattern", "practice"
#   opsin    name OPSIN should read as the same formula; False to skip
#   coords   {atom index: (x, y)} pins heavy atoms (y up, 1 = one bond length)
#   like     id of a resonance form to copy the layout from (same atom order)
#   polar    True/False: what the page says about the molecule's polarity (checked)

S = []


def L(id, formula, title, smiles, note="", role="worked", opsin=None, coords=None, like=None,
      rotate=0, flip=False, polar=None, polar_note=None):
    S.append(dict(id=id, formula=formula, title=title, smiles=smiles, note=note, role=role,
                  opsin=title if opsin is None else opsin, coords=coords, like=like,
                  rotate=rotate, flip=flip, polar=polar, polar_note=polar_note))


# Polarity: set polar=True/False on any structure the text calls polar or nonpolar.
# The build works out polarity from the structure and stops if the two disagree.
# polar_note marks a real exception to the course rules and says why.


LINE3 = {0: (-1, 0), 1: (0, 0), 2: (1, 0)}
PYRAMID = {0: (0, 0), 1: (-0.87, -0.5), 2: (0.87, -0.5), 3: (0, -1)}   # NH3, H3O+: lone pair on top
PLUS = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1)}       # CH4, NH4+
TRIGONAL = {0: (-0.87, -0.5), 1: (0, 0), 2: (0, 1), 3: (0.87, -0.5)}  # O, center, O, O

# ---------------- 2. Drawing Lewis structures ----------------
L("ch2o", "CH2O", "formaldehyde", "C=O", "methanal; preserves tissue samples", coords={0: (0, 0), 1: (0, 1)},
  polar=True)
L("h2o", "H2O", "water", "O", polar=True)
L("nh3", "NH3", "ammonia", "N", coords=PYRAMID, polar=True)
L("co2", "CO2", "carbon dioxide", "O=C=O", polar=False)
L("hcn", "HCN", "hydrogen cyanide", "C#N", coords={0: (0, 0), 1: (1, 0)})
L("c2h4", "C2H4", "ethene", "C=C", "ethylene; ripens fruit", coords={0: (0, 0), 1: (1, 0)})
L("ch4", "CH4", "methane", "C", "a greenhouse gas", role="ref", coords=PLUS, polar=False)
L("n2", "N2", "nitrogen", "N#N", "78% of the air", role="practice", opsin="dinitrogen")
L("ch3oh", "CH3OH", "methanol", "CO", "wood alcohol", role="practice", coords={0: (0, 0), 1: (1, 0)})
L("ch3cl", "CH3Cl", "chloromethane", "CCl", role="practice", coords={0: (0, 0), 1: (1, 0)})
L("h2co3", "H2CO3", "carbonic acid", "OC(O)=O", "forms when CO₂ dissolves in blood", role="practice",
  coords={0: (-0.87, -0.5), 1: (0, 0), 2: (0.87, -0.5), 3: (0, 1)})

# ---------------- 3. Bonding patterns and formal charge ----------------
# Pattern fragments: * is a bond to some other atom. Only the center atom is checked.
# Layouts: atom 1 is the center; bonds point left, up, right, down; lone pairs fill the rest.
X1 = {1: (0, 0), 0: (-1, 0)}
X2 = {1: (0, 0), 0: (-1, 0), 2: (1, 0)}
X3 = {1: (0, 0), 0: (-1, 0), 2: (0, -1), 3: (1, 0)}
X4 = {1: (0, 0), 0: (-1, 0), 2: (0, 1), 3: (1, 0), 4: (0, -1)}
L("pat_c", None, "C, neutral", "*C(*)(*)*", role="pattern", opsin=False, coords=X4)
L("pat_cm", None, "C, −1", "*[C-](*)*", role="pattern", opsin=False, coords=X3)
L("pat_n", None, "N, neutral", "*N(*)*", role="pattern", opsin=False, coords=X3)
L("pat_np", None, "N, +1", "*[N+](*)(*)*", role="pattern", opsin=False, coords=X4)
L("pat_nm", None, "N, −1", "*[N-]*", role="pattern", opsin=False, coords=X2)
L("pat_o", None, "O, neutral", "*O*", role="pattern", opsin=False, coords=X2)
L("pat_op", None, "O, +1", "*[O+](*)*", role="pattern", opsin=False, coords=X3)
L("pat_om", None, "O, −1", "*[O-]", role="pattern", opsin=False, coords=X1)

L("nh4", "NH4^+", "ammonium", "[NH4+]", "N with four bonds: still an octet", coords=PLUS)
L("h3o", "H3O^+", "hydronium", "[OH3+]", "O with three bonds", coords=PYRAMID)
L("oh", "OH^-", "hydroxide", "[OH-]", "O with one bond", coords={0: (0, 0)})
L("co", "CO", "carbon monoxide", "[C-]#[O+]", "binds hemoglobin's iron", role="health", coords={0: (0, 0), 1: (1, 0)})
L("n2o", "N2O", "nitrous oxide", "N#[N+][O-]", "laughing gas; a greenhouse gas", coords=LINE3)
L("cn", "CN^-", "cyanide", "[C-]#N", role="practice", coords={0: (0, 0), 1: (1, 0)})
L("mea", "CH3NH3^+", "methylammonium", "C[NH3+]", role="practice", coords={0: (0, 0), 1: (1, 0)})

# ---------------- 4. Resonance ----------------
L("o3a", "O3", "ozone", "O=[O+][O-]", "form 1", coords={0: (-0.87, -0.5), 1: (0, 0), 2: (0.87, -0.5)},
  polar=True, polar_note="all O–O bonds, but the central O differs from the outer two")
L("o3b", "O3", "ozone", "[O-][O+]=O", "form 2", like="o3a")
L("co3a", "CO3^2-", "carbonate", "[O-]C(=O)[O-]", coords=TRIGONAL)
L("co3b", "CO3^2-", "carbonate", "O=C([O-])[O-]", like="co3a")
L("co3c", "CO3^2-", "carbonate", "[O-]C([O-])=O", like="co3a")
L("hco3a", "HCO3^-", "hydrogen carbonate", "OC(=O)[O-]", "bicarbonate", role="health",
  opsin="hydrogencarbonate", coords=TRIGONAL)
L("hco3b", "HCO3^-", "hydrogen carbonate", "OC([O-])=O", "bicarbonate", role="health",
  opsin="hydrogencarbonate", like="hco3a")
L("ac_a", "CH3COO^-", "acetate", "CC(=O)[O-]", coords=TRIGONAL)
L("ac_b", "CH3COO^-", "acetate", "CC([O-])=O", like="ac_a")
L("no3a", "NO3^-", "nitrate", "[O-][N+](=O)[O-]", role="practice", coords=TRIGONAL)
L("no3b", "NO3^-", "nitrate", "O=[N+]([O-])[O-]", role="practice", like="no3a")
L("no3c", "NO3^-", "nitrate", "[O-][N+]([O-])=O", role="practice", like="no3a")

# ---------------- Page 2: polarity examples ----------------
L("ccl4", "CCl4", "tetrachloromethane", "ClC(Cl)(Cl)Cl", "carbon tetrachloride",
  coords={1: (0, 0), 0: (-1, 0), 2: (0, 1), 3: (1, 0), 4: (0, -1)}, polar=False)
L("chcl3", "CHCl3", "trichloromethane", "ClC(Cl)Cl", "chloroform, an early surgical anesthetic",
  coords={1: (0, 0), 0: (-1, 0), 2: (0, 1), 3: (1, 0)}, polar=True)

# ---------------- Page 2: intermolecular forces ----------------
L("etoh", "C2H6O", "ethanol", "CCO", "boils at 78 °C", polar=True,
  coords=LINE3)
L("dme", "C2H6O", "methoxymethane", "COC", "dimethyl ether; boils at −24 °C", polar=True,
  coords=LINE3)
L("pyr_acid", "C3H4O3", "pyruvic acid", "CC(=O)C(=O)O", "the form in a bottle", opsin="2-oxopropanoic acid",
  coords={0: (-1.73, -0.5), 1: (-0.87, 0), 2: (-0.87, 1), 3: (0, -0.5), 4: (0, -1.5), 5: (0.87, 0)})
L("pyruvate", "C3H3O3^-", "pyruvate", "CC(=O)C(=O)[O-]", "the form in your cells", role="health",
  coords={0: (-1.73, -0.5), 1: (-0.87, 0), 2: (-0.87, 1), 3: (0, -0.5), 4: (0, -1.5), 5: (0.87, 0)})
L("ala", "C3H7NO2", "alanine", "CC([NH3+])C(=O)[O-]", "at body pH: a zwitterion", role="health",
  opsin="2-azaniumylpropanoate",
  coords={0: (-0.87, -0.5), 1: (0, 0), 2: (0, 1), 3: (0.87, -0.5), 4: (0.87, -1.5), 5: (1.73, 0)})

# ---------------- Page 2: beyond the octet ----------------
L("sf6", "SF6", "sulfur hexafluoride", "FS(F)(F)(F)(F)F", "a greenhouse gas", polar=False,
  coords={1: (0, 0), 0: (0, 1), 2: (0.87, 0.5), 3: (0.87, -0.5), 4: (0, -1), 5: (-0.87, -0.5), 6: (-0.87, 0.5)})
L("so4", "SO4^2-", "sulfate", "O=S(=O)([O-])[O-]", coords={1: (0, 0), 0: (0, 1), 2: (1, 0), 3: (0, -1), 4: (-1, 0)})
L("po4", "PO4^3-", "phosphate", "O=P([O-])([O-])[O-]", coords={1: (0, 0), 0: (0, 1), 2: (1, 0), 3: (0, -1), 4: (-1, 0)})
L("xeof4", "XeOF4", "xenon oxide tetrafluoride", "O=[Xe](F)(F)(F)F", opsin="xenon tetrafluoride oxide", polar=True,
  coords={1: (0, 0), 0: (0, 1), 2: (0.94, 0.34), 3: (0.64, -0.77), 4: (-0.64, -0.77), 5: (-0.94, 0.34)})

# Step drawings in "Walk through it" boxes use [[steps id]]: neutral molecules only.

# Practice: (id, structure ids for the answer, explanation)
# The prompt is generated from the first structure's formula and name.
PRACTICE = [
    ("p_n2", ["n2"], "Each N has 5 valence electrons: one pair and three single electrons. One single bond uses one "
     "single electron from each N, leaving two on each. Pair those up as two more bonds: a triple bond. Each N "
     "keeps its pair, so each has 8."),
    ("p_ch3oh", ["ch3oh"], "Connect C and O first, then put the three H on C and the last H on O. Every atom "
     "already has its octet (or duet for H) with single bonds, so the two leftover pairs stay on O."),
    ("p_ch3cl", ["ch3cl"], "C goes in the center because it makes the most bonds. Cl, like F, has one single "
     "electron and three pairs, so it makes one bond and keeps three lone pairs."),
    ("p_h2co3", ["h2co3"], "C is central with three O around it. Two O carry an H. After single bonds, C has one "
     "single electron left and the O without an H has one: pair them into a C=O double bond."),
    ("p_cn", ["cn"], "4 + 5 + 1 for the − charge. A C≡N triple bond leaves one lone pair on each atom. C has three "
     "bonds and a lone pair, one fewer bond than usual, so C is −1. N has its usual three bonds and lone pair, so it "
     "has no charge. Like carbon monoxide, cyanide binds iron through its carbon, which is why it is so toxic."),
    ("p_mea", ["mea"], "4 + 5 + 6 for the six H, minus 1 for the + charge. N bonds to C and three H, with no lone "
     "pair left. Four bonds is one more than N's usual three, so N is +1. Many medicines contain an N like this one "
     "when they're in the body."),
    ("p_no3", ["no3a", "no3b", "no3c"], "5 + 3 × 6, plus 1 for the − charge. N goes in the center with three O "
     "around it. N can have only 8 electrons, so it makes one double bond and two single bonds and keeps no lone "
     "pair. Four bonds make N +1; each singly bonded O has three lone pairs and is −1. Total: +1 − 1 − 1 = −1. "
     "The double bond can go to any of the three O, so nitrate has three resonance forms, just like carbonate."),
]
