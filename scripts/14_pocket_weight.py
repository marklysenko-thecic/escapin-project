"""Phase 4.4: how much weight can each pocket residue carry?

For every residue near the crystal L-Phe position, report:
  - min distance to the L-Phe atoms (after FAD superposition)
  - mean pLDDT of its side-chain atoms (AF3 stores pLDDT in the B-factor column)
  - the superposition error, as context for the distances
"""
import glob
from Bio.PDB import PDBParser, MMCIFParser, Superimposer, NeighborSearch

CUTOFF = 6.0
TARGETS = {
    "escapin":     "models/af3/fold_escapin_acalifornica_fad/*model_0.cif",
    "inktoxin1":   "models/af3/fold_inktoxin1_apunctata_fad/*model_0.cif",
    "aplysianinA": "models/af3/fold_aplysianina_acalifornica_fad/*model_0.cif",
}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
ACIDIC = {"ASP", "GLU"}

xtal = PDBParser(QUIET=True).get_structure("x", "models/reference/2IID.pdb")[0]["A"]
lig = [a for r in xtal if r.get_resname() == "PHE" and r.id[0] != " " for a in r]
xfad = next(r for r in xtal if r.get_resname() == "FAD")

cif = MMCIFParser(QUIET=True)
for name, pat in TARGETS.items():
    st = cif.get_structure(name, glob.glob(pat)[0])[0]
    mfad = next(r for ch in st for r in ch if r.get_resname() == "FAD")
    shared = [a.get_id() for a in xfad if a.get_id() in {b.get_id() for b in mfad}]
    sup = Superimposer()
    sup.set_atoms([xfad[i] for i in shared], [mfad[i] for i in shared])
    sup.apply([a for ch in st for r in ch for a in r])

    prot = [a for ch in st for r in ch if r.id[0] == " " for a in r]
    ns = NeighborSearch(prot)
    hits = {}
    for a in lig:
        for b in ns.search(a.coord, CUTOFF):
            hits.setdefault(b.get_parent(), [])
    for r in hits:
        hits[r] = min(float(((a.coord - l.coord) ** 2).sum()) ** 0.5
                      for a in r for l in lig)

    print(f"\n=== {name}   (FAD superposition error {sup.rms:.2f} A) ===")
    print(f"{'residue':>10} {'dist':>7} {'pLDDT':>7} {'side-chain pLDDT':>18}")
    for r, d in sorted(hits.items(), key=lambda kv: kv[1]):
        allb = [a.get_bfactor() for a in r]
        side = [a.get_bfactor() for a in r if a.get_id() not in BACKBONE] or allb
        flag = "  <-- acidic" if r.get_resname() in ACIDIC else ""
        print(f"{r.get_resname()}{r.id[1]:<6} {d:7.2f} {sum(allb)/len(allb):7.1f} "
              f"{sum(side)/len(side):18.1f}{flag}")
