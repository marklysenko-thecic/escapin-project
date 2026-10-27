"""Phase 4.6b: mean electrostatic potential inside each substrate pocket.

Probe points = the atoms of L-Phe from crystal 2IID, transferred into each
model's frame by superimposing the crystal FAD onto the model FAD.
The potential (kT/e) is read from the APBS .dx map by trilinear interpolation.
"""
import glob, copy
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Superimposer

MODELS = {
    "escapin":     "models/af3/fold_escapin_acalifornica_fad/*model_0.cif",
    "inktoxin1":   "models/af3/fold_inktoxin1_apunctata_fad/*model_0.cif",
    "aplysianinA": "models/af3/fold_aplysianina_acalifornica_fad/*model_0.cif",
    "snake":       "models/af3/fold_laao_crhodostoma_fad/*model_0.cif",
}

def read_dx(path):
    with open(path) as f:
        lines = f.readlines()
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
    grid = np.array(vals).reshape(counts)
    step = np.array([delta[0][0], delta[1][1], delta[2][2]])
    return grid, origin, step

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

xtal = PDBParser(QUIET=True).get_structure("x", "models/reference/2IID.pdb")[0]["A"]
xfad = next(r for r in xtal if r.get_resname() == "FAD")
xphe = [r for r in xtal if r.get_resname() == "PHE" and r.id[0] != " "]
cif = MMCIFParser(QUIET=True)

rows = []
for name, pat in MODELS.items():
    hits = glob.glob(pat)
    if not hits:
        print(f"{name}: model missing"); continue
    m = cif.get_structure(name, hits[0])[0]
    mfad = next(r for ch in m for r in ch if r.get_resname() == "FAD")
    shared = [a.get_id() for a in xfad if a.get_id() in {b.get_id() for b in mfad}]

    # move a COPY of the crystal ligand into this model's frame
    fad_copy = copy.deepcopy(xfad)
    phe = copy.deepcopy(xphe)
    sup = Superimposer()
    sup.set_atoms([mfad[i] for i in shared], [fad_copy[i] for i in shared])
    probes = [a for r in phe for a in r]
    sup.apply(probes)
    pts = np.array([a.coord for a in probes])

    for ph in (5, 8):
        dx = sorted(glob.glob(f"results/apbs/{name}_pH{ph}*.dx"))
        if not dx:
            print(f"{name} pH{ph}: no .dx file"); continue
        grid, origin, step = read_dx(dx[0])
        v = np.array([interp(grid, origin, step, p) for p in pts])
        v = v[~np.isnan(v)]
        rows.append((name, ph, v.mean(), v.min(), v.max(), len(v), sup.rms))

print(f"\n{'protein':12s} {'pH':>3} {'mean kT/e':>10} {'min':>8} {'max':>8} {'n':>4} {'FAD rms':>8}")
for r in rows:
    print(f"{r[0]:12s} {r[1]:>3} {r[2]:10.2f} {r[3]:8.2f} {r[4]:8.2f} {r[5]:>4} {r[6]:8.2f}")
