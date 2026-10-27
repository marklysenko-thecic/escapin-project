"""Phase 7.7: pocket electrostatics for XP_005095149 and its Y207E mutant.

XP_005095149 is the one Aplysia LAAO carrying tyrosine, not glutamate, at the
pocket position that H1 rests on, and it branches off before both the ink-type
and albumen-type groups. It has no measured activity, so this is a prediction
rather than an explanation of known data.

Prediction recorded before running (Phase 7 notes):
  wild type   less negative than the ink-type enzymes (which sit near -18.6)
  Y207E       shifted toward them

The mutant is the matched control. The wild type alone cannot isolate position
207, since it differs from escapin at many positions besides that one.

Runs five AlphaFold models per variant at pH 8, the same pipeline as script 19,
and reports them beside the earlier results.
"""
import csv, copy, glob, shutil, subprocess, sys
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, Superimposer, NeighborSearch

PH = 8
WORK = Path("results/apbs_outlier"); WORK.mkdir(parents=True, exist_ok=True)
CSV = WORK / "outlier_potential.csv"
VARIANTS = {"XP005095149_wt": None, "XP005095149_Y207E": None}
REFERENCE = {"ink toxin 1": -18.6, "aplysianin A": -8.4,
             "escapin": +3.0, "snake model": +5.4, "snake crystal": +8.3}

for t in ("pdb2pqr", "apbs"):
    if shutil.which(t) is None:
        sys.exit(f"'{t}' not found - run: conda activate escapin")

# locate each variant's AF3 folder, whatever case the server used
for name in list(VARIANTS):
    hits = [d for d in glob.glob("models/af3/fold_*")
            if name.lower().replace("_", "") in Path(d).name.lower().replace("_", "")]
    if not hits:
        sys.exit(f"no AF3 folder found for {name} - check models/af3/")
    VARIANTS[name] = hits[0]
    print(f"{name}: {Path(hits[0]).name}")

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
    return (np.array(vals).reshape(counts), origin,
            np.array([delta[0][0], delta[1][1], delta[2][2]]))

def interp(grid, origin, step, p):
    f = (p - origin) / step
    i0 = np.floor(f).astype(int)
    if np.any(i0 < 0) or np.any(i0 + 1 >= np.array(grid.shape)):
        return np.nan
    d = f - i0
    x, y, z = i0
    c = grid[x:x+2, y:y+2, z:z+2]
    c = c[0]*(1-d[0]) + c[1]*d[0]
    c = c[0]*(1-d[1]) + c[1]*d[1]
    return c[0]*(1-d[2]) + c[1]*d[2]

class ProteinOnly(Select):
    def accept_residue(self, r): return r.id[0] == " "

xtal = PDBParser(QUIET=True).get_structure("x", "models/reference/2IID.pdb")[0]["A"]
xfad = next(r for r in xtal if r.get_resname() == "FAD")
xphe = [r for r in xtal if r.get_resname() == "PHE" and r.id[0] != " "]
cif, io = MMCIFParser(QUIET=True), PDBIO()

done = {}
if CSV.exists():
    for row in csv.DictReader(open(CSV)):
        done[(row["variant"], int(row["model"]))] = row

for name, folder in VARIANTS.items():
    for k in range(5):
        if (name, k) in done:
            print(f"{name} m{k}: already done"); continue
        hits = glob.glob(f"{folder}/*model_{k}.cif")
        if not hits:
            print(f"{name} m{k}: cif missing"); continue
        tag = f"{name}_m{k}"
        st = cif.get_structure(tag, hits[0])
        io.set_structure(st); io.save(str(WORK / f"{tag}.pdb"), ProteinOnly())
        r = subprocess.run(["pdb2pqr", "--ff=AMBER", "--titration-state-method=propka",
                            f"--with-ph={PH}", "--apbs-input", f"{tag}.in",
                            f"{tag}.pdb", f"{tag}.pqr"],
                           cwd=WORK, capture_output=True, text=True)
        if not (WORK / f"{tag}.pqr").exists():
            print(f"{tag}: pdb2pqr failed\n{r.stderr[-200:]}"); continue
        print(f"{tag}: running APBS ...", flush=True)
        with open(WORK / f"{tag}.apbs.log", "w") as log:
            subprocess.run(["apbs", f"{tag}.in"], cwd=WORK, stdout=log, stderr=log)
        dx = sorted(WORK.glob(f"{tag}*.dx"))
        if not dx:
            print(f"{tag}: no map written, see {tag}.apbs.log"); continue

        m = st[0]
        mfad = next(rr for ch in m for rr in ch if rr.get_resname() == "FAD")
        shared = [a.get_id() for a in xfad if a.get_id() in {b.get_id() for b in mfad}]
        fc, ph = copy.deepcopy(xfad), copy.deepcopy(xphe)
        sup = Superimposer(); sup.set_atoms([mfad[i] for i in shared], [fc[i] for i in shared])
        atoms = [a for rr in ph for a in rr]; sup.apply(atoms)
        centre = np.mean([a.coord for a in atoms], axis=0)
        ns = NeighborSearch([a for ch in m for rr in ch if rr.id[0] == " " for a in rr])
        g = np.arange(-6, 6.01, 0.8)
        free = [centre + np.array([x, y, z]) for x in g for y in g for z in g
                if x*x + y*y + z*z <= 36 and not ns.search(centre + np.array([x, y, z]), 2.6)]
        grid, origin, step = read_dx(dx[0])
        v = np.array([interp(grid, origin, step, p) for p in free]); v = v[~np.isnan(v)]
        netq = sum(float(l.split()[-2]) for l in open(WORK / f"{tag}.pqr") if l.startswith("ATOM"))
        row = {"variant": name, "model": k, "mean": f"{v.mean():.2f}",
               "sd": f"{v.std():.2f}", "n": len(v), "net_charge": f"{netq:.0f}",
               "fad_rms": f"{sup.rms:.2f}"}
        done[(name, k)] = row
        new = not CSV.exists()
        with open(CSV, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row))
            if new: w.writeheader()
            w.writerow(row)
        print(f"   pocket {v.mean():+.2f} kT/e  (net charge {netq:+.0f})")
        for f in dx: f.unlink()

print(f"\n=== pocket potential at pH {PH} ===")
summary = {}
for name in VARIANTS:
    vals = [float(done[(name, k)]["mean"]) for k in range(5) if (name, k) in done]
    if not vals: continue
    a = np.array(vals); summary[name] = a
    print(f"{name:20s} {' '.join(f'{x:+6.1f}' for x in vals)}  "
          f"mean {a.mean():+6.2f}  sd {a.std(ddof=1):5.2f}")

print("\nfor comparison (from Phase 4):")
for k, v in REFERENCE.items():
    print(f"  {k:16s} {v:+6.1f}")

if len(summary) == 2:
    wt, mu = summary["XP005095149_wt"], summary["XP005095149_Y207E"]
    sep = wt.max() < mu.min() or mu.max() < wt.min()
    print(f"\nY207E shift: {mu.mean() - wt.mean():+.1f} kT/e  "
          f"({'ranges separate' if sep else 'ranges overlap'})")
    print("compare the escapin S206G swap: -13.3 kT/e, ranges separate")
