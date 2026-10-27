"""Phase 4.3: map pocket residues by 3D superposition, not sequence."""
import sys, glob
from Bio.PDB import PDBParser, MMCIFParser, Superimposer, NeighborSearch

SNAKE_POCKET = [89, 90, 208, 223, 227, 322, 372, 374, 430, 464, 465]
TARGETS = {
    "escapin":     "models/af3/fold_escapin_acalifornica_fad/*model_0.cif",
    "inktoxin1":   "models/af3/fold_inktoxin1_apunctata_fad/*model_0.cif",
    "aplysianinA": "models/af3/fold_aplysianina_acalifornica_fad/*model_0.cif",
}

xtal = PDBParser(QUIET=True).get_structure("x", "models/reference/2IID.pdb")[0]["A"]
phe = [r for r in xtal if r.get_resname() == "PHE" and r.id[0] != " "]
lig = [a for r in phe for a in r]

# FAD as the common anchor: superimpose each model's FAD onto the crystal's
xfad = next(r for r in xtal if r.get_resname() == "FAD")

cif = MMCIFParser(QUIET=True)
for name, pat in TARGETS.items():
    st = cif.get_structure(name, glob.glob(pat)[0])[0]
    mfad = next(r for ch in st for r in ch if r.get_resname() == "FAD")
    common = [a.get_id() for a in xfad if a.get_id() in {b.get_id() for b in mfad}]
    sup = Superimposer()
    sup.set_atoms([xfad[i] for i in common], [mfad[i] for i in common])
    sup.apply([a for ch in st for r in ch for a in r])
    print(f"\n{name}: FAD superposition RMSD {sup.rms:.2f} A ({len(common)} atoms)")

    prot = [a for ch in st for r in ch if r.id[0] == " " for a in r]
    ns = NeighborSearch(prot)
    near = {}
    for a in lig:
        for b in ns.search(a.coord, 5.0):
            r = b.get_parent()
            near[r.id[1]] = r.get_resname()
    print(f"  residues within 5 A of the crystal L-Phe position ({len(near)}):")
    print("   " + ", ".join(f"{v}{k}" for k, v in sorted(near.items())))
