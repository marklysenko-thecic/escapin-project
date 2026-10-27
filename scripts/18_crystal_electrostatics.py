"""Phase 4.7: electrostatics of the 2IID crystal structure as a reference.

No modelling and no superposition: the crystal L-Phe sits in the crystal's
own frame, so the pocket location is exact. Compare the result with the
snake AlphaFold model (pH 5: +45.0, pH 8: +4.6 kT/e) to estimate how much
of the model-derived potential is prediction artefact.
"""
import shutil, subprocess, sys
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch

WORK = Path("results/apbs")
for tool in ("pdb2pqr", "apbs"):
    if shutil.which(tool) is None:
        sys.exit(f"'{tool}' not found - run: conda activate escapin")

def read_dx(path):
    lines = open(path).readlines()
    origin, delta, counts, vals, i = None, [], None, [], 0
    while i < len(lines):
        t = lines[i].split()
        if t and t[0] == "object" and "gridpositions" in lines[i]:
            counts = tuple(int(x) for x in t[-3:])
        elif t and t[0] == "origin":
            origin = np.array([float(x) for x in t[1:4]])
        elif t and t[0] == "delta":
            delta.append([float(x) for x in t[1:4]])
        elif t and t[0] == "object" and "array" in lines[i]:
            i += 1
            while i < len(lines) and not lines[i].startswith(("attribute", "object", "component")):
                vals.extend(float(x) for x in lines[i].split())
                i += 1
            continue
        i += 1
    step = np.array([delta[0][0], delta[1][1], delta[2][2]])
    return np.array(vals).reshape(counts), origin, step

def interp(grid, origin, step, p):
    f = (p - origin) / step
    i0 = np.floor(f).astype(int)
    if np.any(i0 < 0) or np.any(i0 + 1 >= np.array(grid.shape)):
        return np.nan
    d = f - i0
    x, y, z = i0
    c = grid[x:x+2, y:y+2, z:z+2]
    c = c[0] * (1 - d[0]) + c[1] * d[0]
    c = c[0] * (1 - d[1]) + c[1] * d[1]
    return c[0] * (1 - d[2]) + c[1] * d[2]

# 1. chain A protein only (no FAD, no ligand, no water)
st = PDBParser(QUIET=True).get_structure("x", "models/reference/2IID.pdb")
chA = st[0]["A"]

class ChainAProtein(Select):
    def accept_chain(self, c):  return c.id == "A"
    def accept_residue(self, r): return r.id[0] == " "

io = PDBIO(); io.set_structure(st)
io.save(str(WORK / "xtal.pdb"), ChainAProtein())
print("wrote results/apbs/xtal.pdb (chain A protein)")

# 2-3. PDB2PQR and APBS, run inside WORK with bare filenames (avoids path bugs)
for ph in (5, 8):
    tag = f"xtal_pH{ph}"
    r = subprocess.run(["pdb2pqr", "--ff=AMBER", "--titration-state-method=propka",
                        f"--with-ph={ph}", "--apbs-input", f"{tag}.in",
                        "xtal.pdb", f"{tag}.pqr"],
                       cwd=WORK, capture_output=True, text=True)
    if not (WORK / f"{tag}.pqr").exists():
        sys.exit(f"pdb2pqr failed at pH {ph}:\n{r.stderr[-400:]}")
    print(f"pH {ph}: pqr written, running APBS (a few minutes)...")
    with open(WORK / f"{tag}.apbs.log", "w") as log:
        subprocess.run(["apbs", f"{tag}.in"], cwd=WORK, stdout=log, stderr=log)
    if not list(WORK.glob(f"{tag}*.dx")):
        sys.exit(f"APBS wrote no map at pH {ph}; see {WORK}/{tag}.apbs.log")

# 4. pocket potential at free-space points around the crystal L-Phe
phe = [r for r in chA if r.get_resname() == "PHE" and r.id[0] != " "]
centre = np.mean([a.coord for r in phe for a in r], axis=0)
prot = [a for r in chA if r.id[0] == " " for a in r]
ns = NeighborSearch(prot)
g = np.arange(-6, 6.01, 0.8)
pts = [centre + np.array([x, y, z]) for x in g for y in g for z in g
       if x*x + y*y + z*z <= 36]
free = np.array([p for p in pts if not ns.search(p, 2.6)])

print(f"\n{'structure':14s} {'pH':>3} {'mean kT/e':>10} {'sd':>7} {'n pts':>6}")
for ph in (5, 8):
    grid, origin, step = read_dx(sorted(WORK.glob(f"xtal_pH{ph}*.dx"))[0])
    v = np.array([interp(grid, origin, step, p) for p in free])
    v = v[~np.isnan(v)]
    print(f"{'2IID crystal':14s} {ph:>3} {v.mean():10.2f} {v.std():7.2f} {len(v):>6}")
print(f"{'snake AF3 model':14s}   5      45.00    5.41    610")
print(f"{'snake AF3 model':14s}   8       4.64    5.31    610")
