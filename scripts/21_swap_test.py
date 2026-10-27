"""Phase 4.10: pocket potential, five models each, including the 206 swap mutants (pH 8).

Purpose: separate model noise from real differences. If each protein's five
models agree with each other but differ between proteins, the gap is real.
If the models of one protein scatter as widely as the gap, it is noise.

Resumable: models already measured are skipped (results kept in the CSV).
"""
import csv, copy, glob, shutil, subprocess, sys
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, Superimposer, NeighborSearch

PH = 8
DELETE_DX = True          # each map is ~60 MB; the numbers are kept in the CSV
WORK = Path("results/apbs_ens"); WORK.mkdir(parents=True, exist_ok=True)
CSV = WORK / "pocket_potential_ensemble.csv"
PROTEINS = {
    "escapin":     "fold_escapin_acalifornica_fad",
    "inktoxin1":   "fold_inktoxin1_apunctata_fad",
    "aplysianinA": "fold_aplysianina_acalifornica_fad",
    "snake":       "fold_laao_crhodostoma_fad",
    "escapin_S206G":   "fold_escapin_s206g_fad",
    "inktoxin1_G206S": "fold_inktoxin1_g206s_fad",
}
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

class ProteinOnly(Select):
    def accept_residue(self, r): return r.id[0] == " "

xtal = PDBParser(QUIET=True).get_structure("x", "models/reference/2IID.pdb")[0]["A"]
xfad = next(r for r in xtal if r.get_resname() == "FAD")
xphe = [r for r in xtal if r.get_resname() == "PHE" and r.id[0] != " "]
cif, io = MMCIFParser(QUIET=True), PDBIO()

done = {}
if CSV.exists():
    for row in csv.DictReader(open(CSV)):
        done[(row["protein"], int(row["model"]))] = row

for name, folder in PROTEINS.items():
    for k in range(5):
        if (name, k) in done:
            print(f"{name} model_{k}: already measured, skipping"); continue
        hits = glob.glob(f"models/af3/{folder}/*model_{k}.cif")
        if not hits:
            print(f"{name} model_{k}: cif not found"); continue
        tag = f"{name}_m{k}"
        st = cif.get_structure(tag, hits[0])
        io.set_structure(st); io.save(str(WORK / f"{tag}.pdb"), ProteinOnly())

        r = subprocess.run(["pdb2pqr", "--ff=AMBER", "--titration-state-method=propka",
                            f"--with-ph={PH}", "--apbs-input", f"{tag}.in",
                            f"{tag}.pdb", f"{tag}.pqr"], cwd=WORK, capture_output=True, text=True)
        if not (WORK / f"{tag}.pqr").exists():
            print(f"{tag}: pdb2pqr failed\n{r.stderr[-300:]}"); continue
        print(f"{tag}: running APBS ...", flush=True)
        with open(WORK / f"{tag}.apbs.log", "w") as log:
            subprocess.run(["apbs", f"{tag}.in"], cwd=WORK, stdout=log, stderr=log)
        dx = sorted(WORK.glob(f"{tag}*.dx"))
        if not dx:
            print(f"{tag}: APBS wrote no map, see {tag}.apbs.log"); continue

        # pocket centre: crystal L-Phe moved into this model's frame via FAD
        m = st[0]
        mfad = next(r for ch in m for r in ch if r.get_resname() == "FAD")
        shared = [a.get_id() for a in xfad if a.get_id() in {b.get_id() for b in mfad}]
        fc, ph = copy.deepcopy(xfad), copy.deepcopy(xphe)
        sup = Superimposer(); sup.set_atoms([mfad[i] for i in shared], [fc[i] for i in shared])
        atoms = [a for rr in ph for a in rr]; sup.apply(atoms)
        centre = np.mean([a.coord for a in atoms], axis=0)

        prot = [a for ch in m for rr in ch if rr.id[0] == " " for a in rr]
        ns = NeighborSearch(prot)
        g = np.arange(-6, 6.01, 0.8)
        pts = [centre + np.array([x, y, z]) for x in g for y in g for z in g
               if x*x + y*y + z*z <= 36]
        free = [p for p in pts if not ns.search(p, 2.6)]

        grid, origin, step = read_dx(dx[0])
        v = np.array([interp(grid, origin, step, p) for p in free]); v = v[~np.isnan(v)]
        netq = sum(float(l.split()[-2]) for l in open(WORK / f"{tag}.pqr") if l.startswith("ATOM"))
        row = {"protein": name, "model": k, "mean": f"{v.mean():.2f}", "sd": f"{v.std():.2f}",
               "n": len(v), "net_charge": f"{netq:.0f}", "fad_rms": f"{sup.rms:.2f}"}
        done[(name, k)] = row
        new = not CSV.exists()
        with open(CSV, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(row)); 
            if new: w.writeheader()
            w.writerow(row)
        print(f"   pocket {v.mean():+.2f} kT/e  (net charge {netq:+.0f})")
        if DELETE_DX:
            for f in dx: f.unlink()

print(f"\n=== pocket potential at pH {PH}, five models each ===")
print(f"{'protein':16s} {'models':>34s} {'mean':>7} {'sd':>6} {'range':>14}")
summary = {}
for name in PROTEINS:
    vals = [float(done[(name, k)]["mean"]) for k in range(5) if (name, k) in done]
    if not vals: continue
    a = np.array(vals); summary[name] = a
    shown = " ".join(f"{x:+6.1f}" for x in vals)
    print(f"{name:16s} {shown:>34s} {a.mean():+7.2f} {a.std(ddof=1) if len(a)>1 else 0:6.2f} "
          f"{a.min():+6.1f}..{a.max():+5.1f}")

def compare(a, b, label):
    if a not in summary or b not in summary:
        print(f"{label}: missing data"); return
    x, y = summary[a], summary[b]
    sep = x.max() < y.min() or y.max() < x.min()
    print(f"{label:34s} {x.mean():+7.1f} -> {y.mean():+7.1f}  shift {y.mean()-x.mean():+6.1f}  "
          f"{'ranges separate' if sep else 'ranges overlap'}")

print("\n=== swap test: does position 206 control the pocket? ===")
compare("escapin", "escapin_S206G", "escapin  S206G (wt -> mutant)")
compare("inktoxin1", "inktoxin1_G206S", "inktoxin1 G206S (wt -> mutant)")
print("\nreference gap escapin vs inktoxin1 (wild types):",
      f"{summary['escapin'].mean() - summary['inktoxin1'].mean():+.1f} kT/e" if 'escapin' in summary and 'inktoxin1' in summary else "n/a")
