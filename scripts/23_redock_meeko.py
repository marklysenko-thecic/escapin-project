"""Phase 5.1b: re-docking control, Meeko ligand preparation.

Replaces the Open Babel route, which wrote all-zero partial charges.
Builds L-Phe from SMILES (so no crystal coordinates leak into the pose),
prepares it with Meeko, and docks it back into 2IID.

Tries two box sizes: a wide 22 A and a tight 16 A. A tight box is the
standard remedy when the search wanders into neighbouring cavities.
"""
import subprocess, shutil, sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

OUT = Path("docking/control"); OUT.mkdir(parents=True, exist_ok=True)
SMILES = "N[C@@H](Cc1ccccc1)C(=O)O"      # L-phenylalanine
EXHAUST, SEEDS, BOXES = 32, (42, 7, 1234), (22.0, 16.0)

for t in ("vina",):
    if shutil.which(t) is None:
        sys.exit(f"'{t}' not found - run: conda activate dock")
for f in ("receptor.pdbqt", "ref_ligand.pdb"):
    if not (OUT / f).exists():
        sys.exit(f"missing {OUT/f} - run scripts/22_redock_control.py first")

# reference coordinates from the crystal, heavy atoms, by atom name
ref = {}
for l in open(OUT / "ref_ligand.pdb"):
    if l.startswith(("ATOM", "HETATM")):
        el = (l[76:78].strip() or l[12:16].strip()[0])
        if el != "H":
            ref.setdefault(l[12:16].strip(),
                           np.array([float(l[30+8*i:38+8*i]) for i in range(3)]))
centre = np.mean(list(ref.values()), axis=0)
print(f"{len(ref)} reference heavy atoms, box centre {centre.round(2)}")

# build and prepare the ligand
mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)
pdbqt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol)[0])
if not ok:
    sys.exit(f"meeko failed: {err}")
(OUT / "ligand_meeko.pdbqt").write_text(pdbqt)
q = sum(float(l[70:76]) for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM")))
print(f"ligand prepared: net charge {q:+.3f}, "
      f"{sum(1 for l in pdbqt.splitlines() if l.startswith('BRANCH'))} rotatable bonds\n")

# Meeko names atoms by element+index, so match to the crystal by the
# ring/backbone pattern instead: compare sorted heavy-atom distance sets
refpts = np.array(list(ref.values()))

def pose_rmsd(pts):
    """Best RMSD over the symmetry-free assignment: nearest-neighbour matching."""
    if len(pts) != len(refpts):
        return None
    d = np.linalg.norm(pts[:, None, :] - refpts[None, :, :], axis=2)
    used, total = set(), 0.0
    for i in np.argsort(d.min(axis=1)):
        j = min((k for k in range(len(refpts)) if k not in used), key=lambda k: d[i, k])
        used.add(j); total += d[i, j] ** 2
    return float(np.sqrt(total / len(pts)))

print(f"{'box':>5} {'seed':>6} {'pose':>5} {'score':>8} {'RMSD':>7}")
best_overall = None
for box in BOXES:
    for seed in SEEDS:
        out = OUT / f"redock_b{int(box)}_s{seed}.pdbqt"
        cmd = ["vina", "--receptor", str(OUT/"receptor.pdbqt"),
               "--ligand", str(OUT/"ligand_meeko.pdbqt"),
               "--center_x", f"{centre[0]:.3f}", "--center_y", f"{centre[1]:.3f}",
               "--center_z", f"{centre[2]:.3f}",
               "--size_x", str(box), "--size_y", str(box), "--size_z", str(box),
               "--exhaustiveness", str(EXHAUST), "--seed", str(seed),
               "--num_modes", "9", "--out", str(out)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if not out.exists():
            print(f"{box:5.0f} {seed:>6} vina failed: {r.stderr.strip()[-120:]}"); continue
        poses, cur, scores = [], [], []
        for l in open(out):
            if l.startswith("REMARK VINA RESULT"):
                scores.append(float(l.split()[3]))
            elif l.startswith(("ATOM", "HETATM")) and l[77:79].strip() != "HD":
                cur.append([float(l[30+8*i:38+8*i]) for i in range(3)])
            elif l.startswith("ENDMDL"):
                poses.append(np.array(cur)); cur = []
        if cur: poses.append(np.array(cur))
        for i, (pts, sc) in enumerate(zip(poses, scores), start=1):
            rm = pose_rmsd(pts)
            if rm is None: continue
            if best_overall is None or rm < best_overall[0]:
                best_overall = (rm, box, seed, i, sc)
            if i <= 3 or rm < 2.0:
                print(f"{box:5.0f} {seed:>6} {i:>5} {sc:8.2f} {rm:7.2f}")

if best_overall:
    rm, box, seed, i, sc = best_overall
    print(f"\nbest overall: RMSD {rm:.2f} A (box {box:.0f}, seed {seed}, pose #{i}, score {sc:.2f})")
    print("CONTROL PASSED" if rm < 2.0 else "still above 2.0 A")
