"""Phase 5.1: re-docking control.

Take L-Phe out of crystal 2IID, dock it back, and measure how far the best
pose lands from where the crystal says it sits. Success = RMSD < 2.0 A.
If this fails, nothing downstream can be trusted, so it runs before any
sea-hare docking.

Outputs:
  docking/control/receptor.pdb     2IID chain A, protein + FAD, no ligand
  docking/control/ref_ligand.pdb   the crystal L-Phe (the answer we hide)
  docking/control/redock_out.pdbqt poses Vina found
"""
import subprocess, sys, shutil
from pathlib import Path
import numpy as np

OUT = Path("docking/control"); OUT.mkdir(parents=True, exist_ok=True)
XTAL = Path("models/reference/2IID.pdb")
EXHAUST, SEED, BOX = 32, 42, 22.0

for tool in ("vina", "obabel"):
    if shutil.which(tool) is None:
        sys.exit(f"'{tool}' not found - run: conda activate dock")
if not XTAL.exists():
    sys.exit(f"missing {XTAL}")

# --- split the crystal: receptor (protein + FAD) vs the ligand we hide ---
rec_lines, lig_lines = [], []
for line in open(XTAL):
    if not line.startswith(("ATOM", "HETATM")):
        continue
    if line[21] != "A":                      # chain A only
        continue
    resn, alt = line[17:20].strip(), line[16]
    if alt not in (" ", "A"):                # drop alternate conformations
        continue
    if resn == "HOH":
        continue
    if resn == "PHE" and line.startswith("HETATM"):
        lig_lines.append(line)
    else:
        rec_lines.append(line)

if not lig_lines:
    sys.exit("no hetero PHE found in chain A - check the crystal file")
(OUT / "receptor.pdb").write_text("".join(rec_lines) + "END\n")
(OUT / "ref_ligand.pdb").write_text("".join(lig_lines) + "END\n")
print(f"receptor {len(rec_lines)} atoms, reference ligand {len(lig_lines)} atoms")

def coords(path, heavy_only=True):
    out = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            el = l[76:78].strip() or l[12:16].strip()[0]
            if heavy_only and el == "H":
                continue
            out.append((l[12:16].strip(), np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    return out

ref = coords(OUT / "ref_ligand.pdb")
centre = np.mean([c for _, c in ref], axis=0)
print(f"box centre {centre.round(2)}, size {BOX} A")

# --- prepare receptor and ligand for Vina ---
subprocess.run(["obabel", str(OUT/"receptor.pdb"), "-O", str(OUT/"receptor.pdbqt"),
                "-xr", "-p", "7.4"], check=True, capture_output=True)
subprocess.run(["obabel", str(OUT/"ref_ligand.pdb"), "-O", str(OUT/"ligand.pdbqt"),
                "-p", "7.4", "--gen3d", "best"], check=True, capture_output=True)
print("prepared receptor.pdbqt and ligand.pdbqt (ligand conformation regenerated)")

# --- dock ---
cmd = ["vina", "--receptor", str(OUT/"receptor.pdbqt"), "--ligand", str(OUT/"ligand.pdbqt"),
       "--center_x", f"{centre[0]:.3f}", "--center_y", f"{centre[1]:.3f}",
       "--center_z", f"{centre[2]:.3f}",
       "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
       "--exhaustiveness", str(EXHAUST), "--seed", str(SEED),
       "--num_modes", "9", "--out", str(OUT/"redock_out.pdbqt")]
print("running vina ...", flush=True)
r = subprocess.run(cmd, capture_output=True, text=True)
(OUT / "vina.log").write_text(r.stdout + r.stderr)
if not (OUT/"redock_out.pdbqt").exists():
    sys.exit(f"vina produced no output; see {OUT}/vina.log\n{r.stderr[-400:]}")

# --- RMSD of each pose to the crystal ligand, matched by atom name ---
refmap = {}
for name, c in ref:
    refmap.setdefault(name, c)

poses, cur, scores = [], [], []
for l in open(OUT/"redock_out.pdbqt"):
    if l.startswith("REMARK VINA RESULT"):
        scores.append(float(l.split()[3]))
    elif l.startswith(("ATOM", "HETATM")):
        el = l[76:78].strip()
        if el != "H":
            cur.append((l[12:16].strip(), np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    elif l.startswith("ENDMDL"):
        poses.append(cur); cur = []
if cur: poses.append(cur)

print(f"\n{'pose':>4} {'score':>8} {'RMSD':>7}  {'matched atoms':>13}")
best = None
for i, (pose, sc) in enumerate(zip(poses, scores), start=1):
    d = [np.linalg.norm(c - refmap[n]) for n, c in pose if n in refmap]
    if not d:
        print(f"{i:>4} {sc:8.2f}    n/a  (atom names did not match)"); continue
    rmsd = float(np.sqrt(np.mean(np.square(d))))
    if best is None or rmsd < best[1]: best = (i, rmsd, sc)
    print(f"{i:>4} {sc:8.2f} {rmsd:7.2f}  {len(d):>13}")

if best:
    i, rmsd, sc = best
    print(f"\nclosest pose: #{i}, RMSD {rmsd:.2f} A, score {sc:.2f} kcal/mol")
    print("CONTROL PASSED" if rmsd < 2.0 else "CONTROL FAILED - fix the setup before docking anything else")
    if i != 1:
        print(f"note: top-scored pose is #1 (RMSD {np.nan if len(poses)<1 else ''}), "
              f"so the scoring function did not rank the correct pose first")
