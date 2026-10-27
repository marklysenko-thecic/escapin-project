"""Phase 7.2: H2 as a CLASSIFICATION question, not a ranking one.

Fig 7A splits cleanly in two: L-Lys (0.42) and L-Arg (0.34) are active, the
other five sit at the no-substrate floor (0.055-0.06). There is no real order
within that floor, so a rank correlation penalises the metric for failing at
something the data cannot define. The answerable question is binary:

    does the salt-bridge metric separate active from inactive substrates?

Tests used, both exact and both appropriate to n = 7:
  Mann-Whitney U  - do the two groups' metric values differ?  (exact, all
                    C(7,2) = 21 splits enumerated)
  AUC             - probability a random active scores above a random inactive;
                    equals U / (n_active * n_inactive), so it is the same
                    statistic read as an effect size
A permutation p-value is computed by enumerating every way to choose which two
of the seven ligands are labelled active.
"""
from itertools import combinations
import glob
import numpy as np

ACTIVITY = {"L-Lys":0.42, "L-Arg":0.34, "D-Lys":0.06, "L-Phe":0.06,
            "L-Met":0.06, "L-Leu":0.055, "L-Trp":0.055}
ACTIVE = {"L-Lys", "L-Arg"}          # the two bars above the floor in Fig 7A
AD = {"A":"C","C":"C","OA":"O","O":"O","N":"N","NA":"N","SA":"S","S":"S"}

def salt_fraction(receptor_pdbqt, run_prefix, ligands):
    ra = []
    for l in open(receptor_pdbqt):
        if l.startswith(("ATOM","HETATM")):
            t = l.split()[-1].upper()
            ra.append((l[17:20].strip(), l[12:16].strip(), AD.get(t, t[0]),
                       np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    n5 = next(c for rn,an,el,c in ra if rn=="FAD" and an=="N5")
    acids = [c for rn,an,el,c in ra if rn in ("ASP","GLU")
             and an in ("OD1","OD2","OE1","OE2") and np.linalg.norm(c-n5) < 14.0]
    out = {}
    for lg in ligands:
        npose = nsalt = 0
        for f in glob.glob(f"docking/runs/{run_prefix}__{lg}__s*.pdbqt"):
            cur = []
            for l in open(f):
                if l.startswith(("ATOM","HETATM")):
                    cur.append(l)
                elif l.startswith("ENDMDL"):
                    ns = [np.array([float(x[30+8*i:38+8*i]) for i in range(3)])
                          for x in cur if AD.get(x.split()[-1].upper(),"X") == "N"]
                    npose += 1
                    if any(min(np.linalg.norm(c-a) for a in acids) < 4.0 for c in ns):
                        nsalt += 1
                    cur = []
        out[lg] = 100*nsalt/npose if npose else float("nan")
    return out

def auc(pos, neg):
    """probability a random positive exceeds a random negative; ties count 0.5"""
    wins = sum((p > n) + 0.5*(p == n) for p in pos for n in neg)
    return wins / (len(pos)*len(neg))

def exact_auc_p(values, labels):
    """P(AUC >= observed) over every way of assigning the same number of actives"""
    k = sum(labels)
    obs = auc([v for v,l in zip(values,labels) if l],
              [v for v,l in zip(values,labels) if not l])
    hits = tot = 0
    for pick in combinations(range(len(values)), k):
        pos = [values[i] for i in pick]
        neg = [values[i] for i in range(len(values)) if i not in pick]
        tot += 1
        if auc(pos, neg) >= obs - 1e-12:
            hits += 1
    return obs, hits/tot, tot

for rec_name, rec_path, prefix in [
        ("ink toxin 1 (A. punctata)", "docking/receptors/inktoxin1.pdbqt", "inktoxin1"),
        ("snake crystal 2IID",        "docking/control/receptor_mkA.pdbqt", "snake_xtal")]:
    ligs = list(ACTIVITY)
    frac = salt_fraction(rec_path, prefix, ligs)
    labels = [1 if l in ACTIVE else 0 for l in ligs]
    vals = [frac[l] for l in ligs]

    print(f"\n=== {rec_name} ===")
    print(f"{'ligand':8s} {'OD':>7} {'active':>7} {'salt %':>8}")
    for l in sorted(ligs, key=lambda z: -frac[z]):
        print(f"{l:8s} {ACTIVITY[l]:7.3f} {'yes' if l in ACTIVE else 'no':>7} {frac[l]:8.0f}")

    a, p, tot = exact_auc_p(vals, labels)
    pos = [v for v,lb in zip(vals,labels) if lb]
    neg = [v for v,lb in zip(vals,labels) if not lb]
    print(f"\n  active:   {sorted(pos, reverse=True)}")
    print(f"  inactive: {sorted(neg, reverse=True)}")
    print(f"  AUC = {a:.3f}   (0.5 = no separation, 1.0 = perfect)")
    print(f"  exact p = {p:.4f}  ({tot} labellings enumerated)")
    print(f"  smallest p reachable with 2 actives of 7: {1/tot:.4f}")
