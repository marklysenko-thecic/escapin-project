"""Phase 7.1: H2 correlation for ink toxin 1 against MEASURED activity.

Activity values are read off Figure 7A (lower panel) of patent WO2004065415A2,
which assays the A. punctata enzyme - the same protein as receptor inktoxin1.
This is the only quantitative substrate data available for any enzyme in the set:
the escapin literature reports only that L-Lys and L-Arg are both effective.

Two correlations are reported because the choice matters:
  (a) all 7 ligands tested in the figure that we also docked
  (b) the same, excluding D-Lys - whose inactivity is stereochemical, something
      the salt-bridge metric provably cannot see (it scored D-Lys like L-Lys)

p-values are exact: every permutation of the activity values is enumerated and
compared against the observed rho, so no distributional assumption is made.
"""
from itertools import permutations
import csv, glob, json, re
from pathlib import Path
import numpy as np

# --- measured activity, OD, patent WO2004065415A2 Fig 7A lower panel ---
# background floor is ~0.06 (the no-substrate level)
ACTIVITY = {
    "L-Lys": 0.42,
    "L-Arg": 0.34,
    "D-Lys": 0.06,
    "L-Phe": 0.06,
    "L-Leu": 0.055,
    "L-Met": 0.06,
    "L-Trp": 0.055,
}
# L-Orn excluded: not tested in the figure, and no escapin data either

AD = {"A":"C","C":"C","OA":"O","O":"O","N":"N","NA":"N","SA":"S","S":"S"}
REC = "docking/receptors/inktoxin1.pdbqt"

def ranks(x):
    x = np.asarray(x, float)
    r = np.empty(len(x))
    order = np.argsort(x)
    r[order] = np.arange(len(x), dtype=float)
    for v in set(x.tolist()):                 # average tied ranks
        m = x == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r

def spearman(a, b):
    ra, rb = ranks(a), ranks(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    den = np.sqrt((ra**2).sum() * (rb**2).sum())
    return float(ra @ rb / den) if den > 0 else 0.0

def exact_p(x, y):
    """one-sided: P(rho >= observed) over all orderings of y"""
    obs = spearman(x, y)
    hits = tot = 0
    for p in permutations(range(len(y))):
        tot += 1
        if spearman(x, [y[i] for i in p]) >= obs - 1e-12:
            hits += 1
    return obs, hits / tot, tot

# --- salt-bridge fraction per ligand, from the stored poses ---
ra = []
for l in open(REC):
    if l.startswith(("ATOM","HETATM")):
        t = l.split()[-1].upper()
        ra.append((l[17:20].strip(), l[12:16].strip(), AD.get(t, t[0]),
                   np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
n5 = next(c for rn, an, el, c in ra if rn == "FAD" and an == "N5")
acids = [c for rn, an, el, c in ra if rn in ("ASP","GLU")
         and an in ("OD1","OD2","OE1","OE2") and np.linalg.norm(c - n5) < 14.0]

frac = {}
for lg in ACTIVITY:
    npose = nsalt = 0
    for f in glob.glob(f"docking/runs/inktoxin1__{lg}__s*.pdbqt"):
        cur = []
        for l in open(f):
            if l.startswith(("ATOM","HETATM")):
                cur.append(l)
            elif l.startswith("ENDMDL"):
                ns = [np.array([float(x[30+8*i:38+8*i]) for i in range(3)])
                      for x in cur if AD.get(x.split()[-1].upper(), "X") == "N"]
                npose += 1
                if any(min(np.linalg.norm(c - a) for a in acids) < 4.0 for c in ns):
                    nsalt += 1
                cur = []
    frac[lg] = 100 * nsalt / npose if npose else float("nan")

ligs = list(ACTIVITY)
print("INK TOXIN 1 (A. punctata) - measured activity vs docking metric\n")
print(f"{'ligand':8s} {'activity OD':>12} {'salt-bridge %':>14}")
for lg in ligs:
    print(f"{lg:8s} {ACTIVITY[lg]:12.3f} {frac[lg]:14.0f}")

for label, keep in [("all 7 ligands", ligs),
                    ("excluding D-Lys", [l for l in ligs if l != "D-Lys"])]:
    x = [frac[l] for l in keep]
    y = [ACTIVITY[l] for l in keep]
    rho, p, tot = exact_p(x, y)
    print(f"\n{label}  (n = {len(keep)})")
    print(f"   Spearman rho = {rho:+.3f}")
    print(f"   exact p      = {p:.4f}   ({tot} permutations enumerated)")
    print(f"   ranks: activity {ranks(y).tolist()}")
    print(f"          metric   {ranks(x).tolist()}")
