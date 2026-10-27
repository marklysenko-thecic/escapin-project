"""Phase 5.2: build the eight docking ligands as zwitterions.

All are L-amino acids except D-Lys, the stereochemistry control.
Side chains are protonated for pH ~8: Lys and Arg charged, His neutral.

  L-Lys, L-Arg   escapin's known substrates (positive side chains)
  L-Phe, L-Leu, L-Met, L-Trp   hydrophobic; preferred by the snake enzyme
  L-Orn          lysine with a shorter side chain - tests chain length
  D-Lys          negative control; LAAOs are L-specific
"""
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

OUT = Path("docking/ligands"); OUT.mkdir(parents=True, exist_ok=True)

LIGANDS = {
    "L-Lys": "[NH3+][C@@H](CCCC[NH3+])C(=O)[O-]",
    "L-Arg": "[NH3+][C@@H](CCCNC(=[NH2+])N)C(=O)[O-]",
    "L-Orn": "[NH3+][C@@H](CCC[NH3+])C(=O)[O-]",
    "D-Lys": "[NH3+][C@H](CCCC[NH3+])C(=O)[O-]",
    "L-Phe": "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]",
    "L-Leu": "[NH3+][C@@H](CC(C)C)C(=O)[O-]",
    "L-Met": "[NH3+][C@@H](CCSC)C(=O)[O-]",
    "L-Trp": "[NH3+][C@@H](Cc1c[nH]c2ccccc12)C(=O)[O-]",
}
EXPECTED_Q = {"L-Lys": 1, "L-Arg": 1, "L-Orn": 1, "D-Lys": 1,
              "L-Phe": 0, "L-Leu": 0, "L-Met": 0, "L-Trp": 0}

print(f"{'ligand':8s} {'formal q':>9} {'net pdbqt':>10} {'heavy':>6} {'rotB':>5} {'MW':>7}")
for name, smi in LIGANDS.items():
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"{name:8s}  SMILES failed"); continue
    fq = Chem.GetFormalCharge(mol)
    mw = Descriptors.MolWt(mol)
    molh = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(molh, randomSeed=42) != 0:
        print(f"{name:8s}  3D embedding failed"); continue
    AllChem.MMFFOptimizeMolecule(molh)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(molh)[0])
    if not ok:
        print(f"{name:8s}  meeko failed: {err}"); continue
    (OUT / f"{name}.pdbqt").write_text(pdbqt)
    q = sum(float(l.split()[-2]) for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM")))
    heavy = sum(1 for l in pdbqt.splitlines()
                if l.startswith(("ATOM", "HETATM")) and l.split()[-1] not in ("HD", "H"))
    rot = sum(1 for l in pdbqt.splitlines() if l.startswith("BRANCH"))
    flag = "" if fq == EXPECTED_Q[name] else "  <-- UNEXPECTED CHARGE"
    print(f"{name:8s} {fq:>+9d} {q:>+10.2f} {heavy:>6} {rot:>5} {mw:7.1f}{flag}")

# stereochemistry check: L and D Lys must be mirror images, not identical
l = Chem.MolFromSmiles(LIGANDS["L-Lys"]); d = Chem.MolFromSmiles(LIGANDS["D-Lys"])
cl = Chem.FindMolChiralCenters(l, includeUnassigned=True)
cd = Chem.FindMolChiralCenters(d, includeUnassigned=True)
print(f"\nchiral centres  L-Lys {cl}   D-Lys {cd}")
print("stereochemistry OK" if cl != cd else "PROBLEM: L and D are the same")
