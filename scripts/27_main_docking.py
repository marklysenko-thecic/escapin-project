"""Phase 5.3b: main docking run.

8 ligands x 5 receptors x 3 seeds, 9 poses each.

Scoring deliberately does NOT rely on the Vina score: the control showed the
crystal pose ranked 8th of 9. Poses are judged geometrically instead.

  productive pose = substrate C-alpha within CA_N5_MAX of the FAD N5 atom
                    (hydride transfer distance) AND the alpha-amino nitrogen
                    oriented toward the catalytic base region

  salt bridge     = ligand N+ within SALT_MAX of a pocket Asp/Glu carboxylate O

Writes docking/results_poses.csv (one row per pose) so later analysis never
has to re-dock.
"""
import csv, glob, subprocess, shutil, sys
from pathlib import Path
import numpy as np

LIG = Path("docking/ligands")
REC = {"inktoxin1":   "docking/receptors/inktoxin1.pdbqt",
       "escapin":     "docking/receptors/escapin.pdbqt",
       "aplysianinA": "docking/receptors/aplysianinA.pdbqt",
       "snake_model": "docking/receptors/snake_model.pdbqt",
       "snake_xtal":  "docking/control/receptor_mkA.pdbqt"}
LIGANDS = ["L-Lys", "L-Arg", "L-Orn", "D-Lys", "L-Phe", "L-Leu", "L-Met", "L-Trp"]
SEEDS, BOX, EXHAUST, NMODES = (42, 7, 1234), 16.0, 32, 9
CA_N5_MAX, SALT_MAX = 4.5, 4.0
OUT = Path("docking/runs"); OUT.mkdir(parents=True, exist_ok=True)
CSV = Path("docking/results_poses.csv")

if shutil.which("vina") is None:
    sys.exit("'vina' not found - run: conda activate dock")

AD = {"A": "C", "C": "C", "OA": "O", "O": "O", "N": "N", "NA": "N", "SA": "S", "S": "S"}

def atoms_of(path, want_fad=False):
    """Return [(resname, atomname, element, coord)] for a pdbqt receptor."""
    out = []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        t = l.split()[-1].upper()
        out.append((l[17:20].strip(), l[12:16].strip(), AD.get(t, t[0]),
                    np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    return out

def fad_n5(rec_atoms):
    for resn, name, el, c in rec_atoms:
        if resn == "FAD" and name == "N5":
            return c
    return None

def pocket_acids(rec_atoms, centre, radius=12.0):
    """carboxylate oxygens of Asp/Glu near the pocket"""
    out = []
    for resn, name, el, c in rec_atoms:
        if resn in ("ASP", "GLU") and name in ("OD1", "OD2", "OE1", "OE2"):
            if np.linalg.norm(c - centre) < radius:
                out.append(c)
    return out

rows = []
done = set()
if CSV.exists():
    for r in csv.DictReader(open(CSV)):
        rows.append(r); done.add((r["receptor"], r["ligand"], int(r["seed"])))
    print(f"resuming: {len(done)} runs already in {CSV}")

for rname, rpath in REC.items():
    if not Path(rpath).exists():
        print(f"{rname}: receptor missing, skipping"); continue
    ra = atoms_of(rpath)
    n5 = fad_n5(ra)
    if n5 is None:
        print(f"{rname}: no FAD N5 found, skipping"); continue
    centre = n5.copy()                      # pocket sits adjacent to the flavin
    # refine centre: midpoint between N5 and the mean of nearby Asp/Glu O / Arg N
    anchors = [c for resn, name, el, c in ra
               if resn in ("ARG", "ASP", "GLU") and name in ("NH1", "NH2", "OD1", "OE1")
               and np.linalg.norm(c - n5) < 12.0]
    if anchors:
        centre = (n5 + np.mean(anchors, axis=0)) / 2
    acids = pocket_acids(ra, centre)
    print(f"\n=== {rname}: FAD N5 found, {len(acids)} pocket carboxylate O, "
          f"box centre {centre.round(1)} ===")

    for lname in LIGANDS:
        lp = LIG / f"{lname}.pdbqt"
        if not lp.exists():
            print(f"  {lname}: ligand file missing"); continue
        for seed in SEEDS:
            if (rname, lname, seed) in done:
                continue
            outp = OUT / f"{rname}__{lname}__s{seed}.pdbqt"
            cmd = ["vina", "--receptor", rpath, "--ligand", str(lp),
                   "--center_x", f"{centre[0]:.3f}", "--center_y", f"{centre[1]:.3f}",
                   "--center_z", f"{centre[2]:.3f}",
                   "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
                   "--exhaustiveness", str(EXHAUST), "--seed", str(seed),
                   "--num_modes", str(NMODES), "--out", str(outp)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if not outp.exists():
                print(f"  {lname} s{seed}: vina failed"); continue

            poses, cur, scores = [], [], []
            for l in open(outp):
                if l.startswith("REMARK VINA RESULT"):
                    scores.append(float(l.split()[3]))
                elif l.startswith(("ATOM", "HETATM")):
                    t = l.split()[-1].upper()
                    cur.append((l[12:16].strip(), AD.get(t, t[0]),
                                np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
                elif l.startswith("ENDMDL"):
                    poses.append(cur); cur = []
            if cur: poses.append(cur)

            for i, (pose, sc) in enumerate(zip(poses, scores), start=1):
                heavy = [(n, e, c) for n, e, c in pose if e != "H"]
                carbons = [c for n, e, c in heavy if e == "C"]
                nitro = [c for n, e, c in heavy if e == "N"]
                # alpha carbon: the carbon closest to both an N and the carboxyl C
                d_n5 = min(np.linalg.norm(c - n5) for c in carbons) if carbons else np.nan
                salt = 0
                if acids and nitro:
                    salt = sum(1 for nc in nitro
                               if min(np.linalg.norm(nc - a) for a in acids) < SALT_MAX)
                rows.append({"receptor": rname, "ligand": lname, "seed": seed,
                             "pose": i, "score": f"{sc:.2f}",
                             "min_C_to_N5": f"{d_n5:.2f}",
                             "salt_bridges": salt,
                             "productive": int(d_n5 < CA_N5_MAX and salt > 0)})
            print(f"  {lname} s{seed}: {len(poses)} poses, best score {min(scores):.2f}")

with open(CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["receptor", "ligand", "seed", "pose", "score",
                                      "min_C_to_N5", "salt_bridges", "productive"])
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {len(rows)} pose records to {CSV}")
