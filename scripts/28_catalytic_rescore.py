"""Phase 5.4: re-score the existing docked poses with a proper catalytic test.

No re-docking: this reads docking/runs/*.pdbqt written by script 27 and
replaces the loose metric with three requirements that must hold together.

Why the old metric failed:
  - "nearest carbon to N5" was 2.8-3.6 A for every ligand, so it never
    discriminated; the productive fraction was really just the salt-bridge count
  - D-Lys scored as high as L-Lys, because distances and elements alone
    cannot see chirality

The catalytic geometry of an L-amino acid oxidase requires:
  1. the ALPHA CARBON (identified by connectivity: bonded to the alpha N and
     to the carboxyl C) within CA_N5_MAX of FAD N5 - this is the hydride donor
  2. the ALPHA AMINO nitrogen within N_BASE_MAX of the catalytic base, so the
     proton can be abstracted
  3. the C-alpha hydride vector pointing at N5 (angle below ANGLE_MAX), which
     is the step that distinguishes L from D

Alpha atoms are found from the ligand SMILES via RDKit, matched onto the pose
by atom order (Meeko preserves it), so no PDB atom names are needed.
"""
import csv, glob, re, sys
from pathlib import Path
import numpy as np
from rdkit import Chem

CA_N5_MAX, N_BASE_MAX, ANGLE_MAX = 4.5, 6.0, 60.0
AD = {"A": "C", "C": "C", "OA": "O", "O": "O", "N": "N", "NA": "N", "SA": "S", "S": "S"}

SMILES = {
    "L-Lys": "[NH3+][C@@H](CCCC[NH3+])C(=O)[O-]",
    "L-Arg": "[NH3+][C@@H](CCCNC(=[NH2+])N)C(=O)[O-]",
    "L-Orn": "[NH3+][C@@H](CCC[NH3+])C(=O)[O-]",
    "D-Lys": "[NH3+][C@H](CCCC[NH3+])C(=O)[O-]",
    "L-Phe": "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]",
    "L-Leu": "[NH3+][C@@H](CC(C)C)C(=O)[O-]",
    "L-Met": "[NH3+][C@@H](CCSC)C(=O)[O-]",
    "L-Trp": "[NH3+][C@@H](Cc1c[nH]c2ccccc12)C(=O)[O-]",
}
# backbone pattern: alpha-N, alpha-C, carboxyl C, the two O
BACKBONE = Chem.MolFromSmarts("[NX4+,NX3][CX4H1]C(=O)[O-,OX2H1]")

def backbone_indices(smi):
    """heavy-atom indices of (alpha N, alpha C) in Meeko's output order"""
    mol = Chem.MolFromSmiles(smi)
    m = mol.GetSubstructMatch(BACKBONE)
    if not m:
        return None
    # Meeko writes heavy atoms in RDKit order after AddHs; heavy indices are
    # unchanged by AddHs, so the SMILES indices map directly onto the heavy list
    return m[0], m[1]

def heavy_from_pose(lines):
    out = []
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            t = l.split()[-1].upper()
            if t in ("HD", "H"):
                continue
            out.append((AD.get(t, t[0]),
                        np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    return out

def receptor_atoms(path):
    out = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            t = l.split()[-1].upper()
            out.append((l[17:20].strip(), l[12:16].strip(), AD.get(t, t[0]), int(l[22:26]),
                        np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    return out

REC_PATH = {"inktoxin1": "docking/receptors/inktoxin1.pdbqt",
            "escapin": "docking/receptors/escapin.pdbqt",
            "aplysianinA": "docking/receptors/aplysianinA.pdbqt",
            "snake_model": "docking/receptors/snake_model.pdbqt",
            "snake_xtal": "docking/control/receptor_mkA.pdbqt"}

cache = {}
def rec_info(name):
    if name in cache:
        return cache[name]
    ra = receptor_atoms(REC_PATH[name])
    n5 = next((c for rn, an, el, ri, c in ra if rn == "FAD" and an == "N5"), None)
    if n5 is None:
        cache[name] = None; return None
    # catalytic base candidates: His/Lys/Arg nitrogens within 10 A of N5
    base = [c for rn, an, el, ri, c in ra
            if el == "N" and rn in ("HIS", "LYS", "ARG") and np.linalg.norm(c - n5) < 10.0]
    acids = [c for rn, an, el, ri, c in ra
             if rn in ("ASP", "GLU") and an in ("OD1", "OD2", "OE1", "OE2")
             and np.linalg.norm(c - n5) < 14.0]
    cache[name] = (n5, base, acids)
    return cache[name]

idx = {k: backbone_indices(v) for k, v in SMILES.items()}
for k, v in idx.items():
    if v is None:
        sys.exit(f"backbone not matched for {k}")

rows = []
for f in sorted(glob.glob("docking/runs/*.pdbqt")):
    m = re.match(r"(.+?)__(.+?)__s(\d+)\.pdbqt", Path(f).name)
    if not m:
        continue
    rname, lname, seed = m.group(1), m.group(2), int(m.group(3))
    info = rec_info(rname)
    if info is None or lname not in idx:
        continue
    n5, bases, acids = info
    iN, iCA = idx[lname]

    poses, cur, scores = [], [], []
    for l in open(f):
        if l.startswith("REMARK VINA RESULT"):
            scores.append(float(l.split()[3]))
        elif l.startswith(("ATOM", "HETATM")):
            cur.append(l)
        elif l.startswith("ENDMDL"):
            poses.append(cur); cur = []
    if cur: poses.append(cur)

    for i, (pl, sc) in enumerate(zip(poses, scores), start=1):
        h = heavy_from_pose(pl)
        if max(iN, iCA) >= len(h):
            continue
        aN, aCA = h[iN][1], h[iCA][1]
        d_ca = float(np.linalg.norm(aCA - n5))
        d_base = min((float(np.linalg.norm(aN - b)) for b in bases), default=np.inf)
        # hydride vector: opposite the mean direction of CA's heavy neighbours
        others = np.array([c for j, (e, c) in enumerate(h) if j != iCA])
        near = others[np.argsort(np.linalg.norm(others - aCA, axis=1))[:3]]
        hvec = aCA - near.mean(axis=0)
        tvec = n5 - aCA
        cosang = float(hvec @ tvec / (np.linalg.norm(hvec) * np.linalg.norm(tvec) + 1e-9))
        ang = float(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
        salt = sum(1 for e, c in h if e == "N"
                   and min((np.linalg.norm(c - a) for a in acids), default=99) < 4.0)
        ok = int(d_ca < CA_N5_MAX and d_base < N_BASE_MAX and ang < ANGLE_MAX)
        rows.append({"receptor": rname, "ligand": lname, "seed": seed, "pose": i,
                     "score": f"{sc:.2f}", "CA_to_N5": f"{d_ca:.2f}",
                     "N_to_base": f"{d_base:.2f}", "hydride_angle": f"{ang:.0f}",
                     "salt_bridges": salt, "catalytic": ok})

out = Path("docking/results_catalytic.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print(f"re-scored {len(rows)} poses -> {out}\n")

LIG = ["L-Lys", "L-Arg", "L-Orn", "D-Lys", "L-Phe", "L-Leu", "L-Met", "L-Trp"]
RECS = ["inktoxin1", "escapin", "aplysianinA", "snake_model", "snake_xtal"]
for title, key in [("CATALYTIC POSE FRACTION (%)", "catalytic")]:
    print(title)
    print(f"{'ligand':8s}" + "".join(f"{r[:11]:>13s}" for r in RECS))
    for lg in LIG:
        line = f"{lg:8s}"
        for rc in RECS:
            s = [r for r in rows if r["ligand"] == lg and r["receptor"] == rc]
            line += f"{100*sum(r[key] for r in s)/len(s):13.0f}" if s else f"{'-':>13}"
        print(line)

print("\nSPREAD OF CA-to-N5 (A): min / median / max, all receptors pooled")
for lg in LIG:
    v = sorted(float(r["CA_to_N5"]) for r in rows if r["ligand"] == lg)
    print(f"  {lg:8s} {v[0]:5.2f} {v[len(v)//2]:6.2f} {v[-1]:6.2f}")

print("\nL-Lys vs D-Lys, catalytic fraction by receptor:")
for rc in RECS:
    a = [r for r in rows if r["receptor"] == rc and r["ligand"] == "L-Lys"]
    b = [r for r in rows if r["receptor"] == rc and r["ligand"] == "D-Lys"]
    if a and b:
        fa = 100*sum(r["catalytic"] for r in a)/len(a)
        fb = 100*sum(r["catalytic"] for r in b)/len(b)
        print(f"  {rc:12s} L {fa:5.0f}%   D {fb:5.0f}%   difference {fa-fb:+5.0f}")
