"""Phase 5.6: final docking analysis on validated metrics only.

Metrics kept:
  salt_fraction  - % of poses where a ligand nitrogen sits within 4 A of a
                   pocket Asp/Glu carboxylate oxygen. Validated indirectly:
                   it tracks side-chain charge as H1 predicts.
  mean_salt      - mean number of such contacts per pose
  min_CA_to_N5   - closest approach of the true alpha carbon to FAD N5,
                   using the coordinate-verified atom indices

Metrics dropped, with reasons recorded in the output:
  vina score         - ranked the crystal pose 8th of 9 in the control
  catalytic geometry - requires the alpha-amino N within ~6 A of the catalytic
                       base, which the 2IID crystal itself cannot satisfy
                       (His223 NE2 is 9.8 A from FAD N5)

Correlates the salt metric against the published activity ordering and reports
Spearman rho with an exact permutation p-value.
"""
import csv, glob, json, re
from itertools import permutations
from pathlib import Path
import numpy as np

AD = {"A":"C","C":"C","OA":"O","O":"O","N":"N","NA":"N","SA":"S","S":"S"}
REC = {"inktoxin1":"docking/receptors/inktoxin1.pdbqt",
       "escapin":"docking/receptors/escapin.pdbqt",
       "aplysianinA":"docking/receptors/aplysianinA.pdbqt",
       "snake_model":"docking/receptors/snake_model.pdbqt",
       "snake_xtal":"docking/control/receptor_mkA.pdbqt"}
IDX = json.load(open("docking/ligand_atom_index.json"))
LIG = ["L-Lys","L-Arg","L-Orn","D-Lys","L-Phe","L-Leu","L-Met","L-Trp"]

# published activity, higher = more active.
# escapin (Yang 2005) and A. punctata APIT (patent WO2004065415 Fig 7A):
# high on L-Lys and L-Arg, background on the hydrophobic set, D-Lys inactive.
ACT_SEAHARE = {"L-Lys":4,"L-Arg":3,"L-Orn":2,"L-Phe":1,"L-Leu":1,"L-Met":1,"L-Trp":1,"D-Lys":0}
# snake LAAO prefers bulky hydrophobics (Phe, Trp, Leu, Met) over basic residues
ACT_SNAKE   = {"L-Phe":4,"L-Trp":3,"L-Leu":2,"L-Met":2,"L-Lys":1,"L-Arg":1,"L-Orn":1,"D-Lys":0}

def rank(x):
    o = np.argsort(np.argsort(x)).astype(float)
    # average ties
    x = np.asarray(x, float)
    for v in set(x):
        m = x == v
        if m.sum() > 1: o[m] = o[m].mean()
    return o

def spearman(a, b):
    ra, rb = rank(a), rank(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    return float(ra @ rb / (np.sqrt((ra**2).sum() * (rb**2).sum()) + 1e-12))

def perm_p(a, b):
    obs = spearman(a, b)
    n = 0; tot = 0
    for p in permutations(range(len(b))):
        tot += 1
        if spearman(a, [b[i] for i in p]) >= obs: n += 1
    return obs, n / tot

def rec_info(p):
    ra = []
    for l in open(p):
        if l.startswith(("ATOM","HETATM")):
            t = l.split()[-1].upper()
            ra.append((l[17:20].strip(), l[12:16].strip(), AD.get(t,t[0]),
                       np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    n5 = next((c for rn,an,el,c in ra if rn=="FAD" and an=="N5"), None)
    acid = [c for rn,an,el,c in ra if rn in ("ASP","GLU")
            and an in ("OD1","OD2","OE1","OE2") and np.linalg.norm(c-n5) < 14.0]
    return n5, acid

cache = {k: rec_info(v) for k, v in REC.items()}
data = {}
for f in sorted(glob.glob("docking/runs/*.pdbqt")):
    m = re.match(r"(.+?)__(.+?)__s(\d+)\.pdbqt", Path(f).name)
    if not m: continue
    rn, lg = m.group(1), m.group(2)
    if rn not in REC or lg not in IDX: continue
    n5, acids = cache[rn]
    iCA = IDX[lg]["CA"]
    poses, cur = [], []
    for l in open(f):
        if l.startswith(("ATOM","HETATM")): cur.append(l)
        elif l.startswith("ENDMDL"): poses.append(cur); cur = []
    if cur: poses.append(cur)
    for pl in poses:
        h = []
        for l in pl:
            t = l.split()[-1].upper()
            if t in ("HD","H"): continue
            h.append((AD.get(t,t[0]), np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
        if iCA >= len(h): continue
        salt = sum(1 for e,c in h if e=="N"
                   and min((np.linalg.norm(c-a) for a in acids), default=99) < 4.0)
        d = float(np.linalg.norm(h[iCA][1] - n5))
        data.setdefault((rn,lg), []).append((salt, d))

print("SALT-BRIDGE FRACTION (% of poses with >=1 contact)")
print(f"{'ligand':8s}" + "".join(f"{r[:11]:>13s}" for r in REC))
frac = {}
for lg in LIG:
    line = f"{lg:8s}"
    for rn in REC:
        v = data.get((rn,lg))
        if v:
            fr = 100*sum(1 for s,_ in v if s>0)/len(v)
            frac[(rn,lg)] = fr
            line += f"{fr:13.0f}"
        else: line += f"{'-':>13}"
    print(line)

print("\nMEAN SALT CONTACTS PER POSE")
print(f"{'ligand':8s}" + "".join(f"{r[:11]:>13s}" for r in REC))
for lg in LIG:
    line = f"{lg:8s}"
    for rn in REC:
        v = data.get((rn,lg))
        line += f"{np.mean([s for s,_ in v]):13.2f}" if v else f"{'-':>13}"
    print(line)

print("\nCLOSEST ALPHA-CARBON APPROACH TO FAD N5 (A)")
print(f"{'ligand':8s}" + "".join(f"{r[:11]:>13s}" for r in REC))
for lg in LIG:
    line = f"{lg:8s}"
    for rn in REC:
        v = data.get((rn,lg))
        line += f"{min(d for _,d in v):13.2f}" if v else f"{'-':>13}"
    print(line)

print("\n=== H2: does the salt metric reproduce the published ordering? ===")
print(f"{'receptor':14s} {'rho':>7} {'p (exact)':>11}  reference")
for rn in REC:
    act = ACT_SNAKE if rn.startswith("snake") else ACT_SEAHARE
    ligs = [l for l in LIG if (rn,l) in frac]
    x = [frac[(rn,l)] for l in ligs]
    y = [act[l] for l in ligs]
    rho, p = perm_p(x, y)
    src = "snake hydrophobic preference" if rn.startswith("snake") else "Yang 2005 / APIT patent"
    print(f"{rn:14s} {rho:+7.2f} {p:11.4f}  {src}")

print("\nSTEREOCHEMISTRY CONTROL (D-Lys should score below L-Lys)")
for rn in REC:
    a, b = frac.get((rn,"L-Lys")), frac.get((rn,"D-Lys"))
    if a is not None and b is not None:
        verdict = "as expected" if a > b + 5 else ("not separated" if abs(a-b) <= 5 else "INVERTED")
        print(f"  {rn:14s} L {a:5.0f}%  D {b:5.0f}%   {verdict}")
print("\nNote: the metric is charge-based and has no stereochemical term,")
print("so failure to separate L from D is expected, not a surprise.")
