"""Phase 4.1: identify pocket residues in 2IID from the bound L-Phe."""
from Bio.PDB import PDBParser, NeighborSearch

p = PDBParser(QUIET=True)
st = p.get_structure("x", "models/reference/2IID.pdb")
ch = st[0]["A"]

# the bound substrate and the cofactor in chain A
phe = [r for r in ch if r.get_resname() == "PHE" and r.id[0] != " "]
fad = [r for r in ch if r.get_resname() == "FAD"]
print(f"ligands in chain A: {len(phe)} PHE (hetero), {len(fad)} FAD")
if not phe:
    print("no hetero PHE found; listing all hetero residues:")
    print(sorted({r.get_resname() for r in ch if r.id[0] != " "}))
    raise SystemExit

lig_atoms = [a for r in phe for a in r]
prot_atoms = [a for r in ch if r.id[0] == " " for a in r]
ns = NeighborSearch(prot_atoms)

near = set()
for a in lig_atoms:
    for b in ns.search(a.coord, 5.0):
        r = b.get_parent()
        near.add((r.id[1], r.get_resname()))

print(f"\npocket residues within 5 A of L-Phe ({len(near)}):")
for num, name in sorted(near):
    print(f"   {name}{num}")
