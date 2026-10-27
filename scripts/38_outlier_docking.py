"""Phase 7.8: docking for XP_005095149 and its Y207E mutant.

Completes the outlier prediction with the second, independent method.
The electrostatics put the wild-type pocket at +13.6 kT/e, the most positive
in the project. If that is real, basic substrates should form fewer salt
bridges here than in any other Aplysia protein.

Prediction, recorded before running:
  L-Lys and L-Arg salt-bridge fractions well below the ink-type values
  (ink toxin 1: L-Lys 100%, L-Arg 100%) and at or below the snake crystal
  (L-Lys 19%, L-Arg 50%).

Same settings as script 27 so the numbers are comparable: 16 A box centred
on the FAD N5 / pocket-anchor midpoint, exhaustiveness 32, three seeds,
nine poses, salt bridge = ligand N within 4 A of a pocket Asp/Glu carboxylate.
"""
import csv, glob, subprocess, shutil, sys
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select

AD = {"A":"C","C":"C","OA":"O","O":"O","N":"N","NA":"N","SA":"S","S":"S"}
LIGANDS = ["L-Lys", "L-Arg", "L-Orn", "D-Lys", "L-Phe", "L-Leu", "L-Met", "L-Trp"]
SEEDS, BOX, EXHAUST, NMODES = (42, 7, 1234), 16.0, 32, 9
VARIANTS = {"XP005095149_wt": None, "XP005095149_Y207E": None}
REC = Path("docking/receptors"); RUNS = Path("docking/runs")
REC.mkdir(parents=True, exist_ok=True); RUNS.mkdir(parents=True, exist_ok=True)
REFERENCE = {"ink toxin 1": (100, 100), "escapin": (85, 100),
             "aplysianin A": (81, 100), "snake crystal": (19, 50)}

for t in ("vina",):
    if shutil.which(t) is None:
        sys.exit(f"'{t}' not found - run: conda activate dock")
MK = shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")
if MK is None:
    sys.exit("meeko receptor tool not found - run: conda activate dock")

for name in list(VARIANTS):
    key = name.lower().replace("_", "")
    hits = [d for d in glob.glob("models/af3/fold_*")
            if key in Path(d).name.lower().replace("_", "")]
    if not hits:
        sys.exit(f"no AF3 folder for {name}")
    VARIANTS[name] = hits[0]
    print(f"{name}: {Path(hits[0]).name}")

class Keep(Select):
    def accept_residue(self, r): return r.id[0] == " " or r.get_resname() == "FAD"

cif, io = MMCIFParser(QUIET=True), PDBIO()
for name, folder in VARIANTS.items():
    pq = REC / f"{name}.pdbqt"
    if pq.exists():
        print(f"{name}: receptor already prepared"); continue
    st = cif.get_structure(name, glob.glob(f"{folder}/*model_0.cif")[0])
    pdb = REC / f"{name}.pdb"
    io.set_structure(st); io.save(str(pdb), Keep())
    r = subprocess.run([MK, "--read_pdb", str(pdb), "-o", str(REC / name), "-p",
                        "--default_altloc", "A"], capture_output=True, text=True)
    if not pq.exists():
        print("   FAILED:\n   " + "\n   ".join((r.stdout + r.stderr).strip().splitlines()[-5:]))
        continue
    q = sum(float(l.split()[-2]) for l in open(pq) if l.startswith(("ATOM", "HETATM")))
    print(f"   {pq.name}: net charge {q:+.2f}")

def rec_geom(path):
    ra = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            t = l.split()[-1].upper()
            ra.append((l[17:20].strip(), l[12:16].strip(), AD.get(t, t[0]),
                       np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    n5 = next((c for rn, an, el, c in ra if rn == "FAD" and an == "N5"), None)
    anchors = [c for rn, an, el, c in ra
               if rn in ("ARG", "ASP", "GLU") and an in ("NH1", "NH2", "OD1", "OE1")
               and np.linalg.norm(c - n5) < 12.0]
    centre = (n5 + np.mean(anchors, axis=0)) / 2 if anchors else n5
    acids = [c for rn, an, el, c in ra if rn in ("ASP", "GLU")
             and an in ("OD1", "OD2", "OE1", "OE2") and np.linalg.norm(c - n5) < 14.0]
    return n5, centre, acids

rows = []
for name in VARIANTS:
    pq = REC / f"{name}.pdbqt"
    if not pq.exists(): continue
    n5, centre, acids = rec_geom(pq)
    print(f"\n=== {name}: {len(acids)} pocket carboxylate O, centre {centre.round(1)} ===")
    for lg in LIGANDS:
        lp = Path(f"docking/ligands/{lg}.pdbqt")
        if not lp.exists():
            print(f"  {lg}: ligand missing"); continue
        npose = nsalt = 0
        for seed in SEEDS:
            out = RUNS / f"{name}__{lg}__s{seed}.pdbqt"
            if not out.exists():
                subprocess.run(["vina", "--receptor", str(pq), "--ligand", str(lp),
                                "--center_x", f"{centre[0]:.3f}", "--center_y", f"{centre[1]:.3f}",
                                "--center_z", f"{centre[2]:.3f}",
                                "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
                                "--exhaustiveness", str(EXHAUST), "--seed", str(seed),
                                "--num_modes", str(NMODES), "--out", str(out)],
                               capture_output=True, text=True)
            if not out.exists():
                continue
            cur = []
            for l in open(out):
                if l.startswith(("ATOM", "HETATM")):
                    cur.append(l)
                elif l.startswith("ENDMDL"):
                    ns = [np.array([float(x[30+8*i:38+8*i]) for i in range(3)])
                          for x in cur if AD.get(x.split()[-1].upper(), "X") == "N"]
                    npose += 1
                    if any(min(np.linalg.norm(c - a) for a in acids) < 4.0 for c in ns):
                        nsalt += 1
                    cur = []
        frac = 100 * nsalt / npose if npose else float("nan")
        rows.append({"receptor": name, "ligand": lg, "poses": npose,
                     "salt_fraction": f"{frac:.0f}"})
        print(f"  {lg:7s} {frac:5.0f}%  ({npose} poses)")

with open("docking/outlier_salt.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["receptor", "ligand", "poses", "salt_fraction"])
    w.writeheader(); w.writerows(rows)

print("\n=== basic substrates, salt-bridge fraction ===")
print(f"{'receptor':22s} {'L-Lys':>7} {'L-Arg':>7}")
for name in VARIANTS:
    g = {r["ligand"]: r["salt_fraction"] for r in rows if r["receptor"] == name}
    if g: print(f"{name:22s} {g.get('L-Lys','-'):>6}% {g.get('L-Arg','-'):>6}%")
for k, (a, b) in REFERENCE.items():
    print(f"{k+' (earlier)':22s} {a:>6}% {b:>6}%")
