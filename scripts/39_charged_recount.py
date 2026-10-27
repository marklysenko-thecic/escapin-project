"""Phase 5.7: recount salt bridges using FORMAL CHARGES, not nitrogen atoms.

Why this revision: the original metric counted any ligand nitrogen within 4 A
of a pocket Asp/Glu carboxylate. L-Trp scored 85-89% in the sea-hare receptors
despite being uncharged, because its indole nitrogen was counted. That N is a
neutral hydrogen-bond donor, not a charge, so those contacts are not salt
bridges.

This script recounts the SAME docked poses - no re-docking - keeping only
nitrogens that carry a positive formal charge:
    L-Lys  side-chain NH3+      (plus the alpha-amino group)
    L-Arg  guanidinium NH2+     (plus the alpha-amino group)
    L-Orn  side-chain NH3+      (plus the alpha-amino group)
    L-Trp  alpha-amino only; the indole N is now excluded
    L-Phe, L-Leu, L-Met  alpha-amino only

The metric was revised after seeing the results. The reason is principled
rather than cosmetic, but the order should be recorded: see the commit message
and the Phase 5 notes.

Charged-atom indices are found by matching RDKit coordinates to the Meeko
output, the same verified approach used in scripts 29 and 30.

Writes docking/results_charged.csv and prints old vs new side by side.
"""
import csv, glob, re, sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

SMILES = {
    "L-Lys": "[NH3+][C@@H](CCCC[NH3+])C(=O)[O-]",
    "L-Arg": "[NH3+][C@@H](CCCNC(=[NH2+])N)C(=O)[O-]",
    "L-Orn": "[NH3+][C@@H](CCC[NH3+])C(=O)[O-]",
    "D-Lys": "[NH3+][C@H](CCCC[NH3+])C(=O)[O-]",
    "L-Phe": "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]",
    "L-Leu": "[NH3+][C@@H](CC(C)C)C(=O)[O-]",
    "L-Met": "[NH3+][C@@H](CCSC)C(=O)[O-]",
    "L-Trp": "[NH3+][C@@H](Cc1c[nH]c2ccccc12)C(=O)[O-]",
}
REC = {"inktoxin1":"docking/receptors/inktoxin1.pdbqt",
       "escapin":"docking/receptors/escapin.pdbqt",
       "aplysianinA":"docking/receptors/aplysianinA.pdbqt",
       "snake_model":"docking/receptors/snake_model.pdbqt",
       "snake_xtal":"docking/control/receptor_mkA.pdbqt",
       "XP005095149_wt":"docking/receptors/XP005095149_wt.pdbqt",
       "XP005095149_Y207E":"docking/receptors/XP005095149_Y207E.pdbqt"}
AD = {"A":"C","C":"C","OA":"O","O":"O","N":"N","NA":"N","SA":"S","S":"S"}
CUTOFF = 4.0

def charged_n_indices(smi):
    """heavy-atom positions (in Meeko output order) of positively charged N"""
    mol = Chem.MolFromSmiles(smi)
    want = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "N" and a.GetFormalCharge() > 0]
    molh = Chem.AddHs(mol)
    AllChem.EmbedMolecule(molh, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(molh)
    conf = molh.GetConformer()
    coords = {i: np.array([conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                           conf.GetAtomPosition(i).z]) for i in want}
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(molh)[0])
    if not ok:
        sys.exit(f"meeko failed for {smi}: {err}")
    heavy = [np.array([float(l[30+8*i:38+8*i]) for i in range(3)])
             for l in pdbqt.splitlines()
             if l.startswith(("ATOM","HETATM")) and l.split()[-1] not in ("HD","H")]
    idx = []
    for i, c in coords.items():
        d = [np.linalg.norm(h - c) for h in heavy]
        j = int(np.argmin(d))
        if d[j] > 0.05:
            sys.exit(f"could not match atom {i} for {smi}")
        idx.append(j)
    return sorted(idx), len(want)

CH = {}
print("charged nitrogens per ligand")
for name, smi in SMILES.items():
    idx, n = charged_n_indices(smi)
    CH[name] = idx
    print(f"  {name:7s} {n} charged N  -> heavy-atom positions {idx}")

def acids_of(path):
    out, n5 = [], None
    for l in open(path):
        if not l.startswith(("ATOM","HETATM")):
            continue
        rn, an = l[17:20].strip(), l[12:16].strip()
        c = np.array([float(l[30+8*i:38+8*i]) for i in range(3)])
        if rn == "FAD" and an == "N5":
            n5 = c
        if rn in ("ASP","GLU") and an in ("OD1","OD2","OE1","OE2"):
            out.append((an, c))
    if n5 is None:
        return []
    return [c for an, c in out if np.linalg.norm(c - n5) < 14.0]

rows = []
for rname, rpath in REC.items():
    if not Path(rpath).exists():
        print(f"{rname}: receptor file missing, skipping"); continue
    acids = acids_of(rpath)
    for lname in SMILES:
        files = sorted(glob.glob(f"docking/runs/{rname}__{lname}__s*.pdbqt"))
        if not files:
            continue
        npose = old_hit = new_hit = 0
        for f in files:
            cur = []
            for l in open(f):
                if l.startswith(("ATOM","HETATM")):
                    cur.append(l)
                elif l.startswith("ENDMDL"):
                    heavy = [(AD.get(x.split()[-1].upper(), "X"),
                              np.array([float(x[30+8*i:38+8*i]) for i in range(3)]))
                             for x in cur if x.split()[-1].upper() not in ("HD","H")]
                    allN = [c for e, c in heavy if e == "N"]
                    chgN = [heavy[i][1] for i in CH[lname] if i < len(heavy)]
                    npose += 1
                    if any(min(np.linalg.norm(c-a) for a in acids) < CUTOFF for c in allN):
                        old_hit += 1
                    if any(min(np.linalg.norm(c-a) for a in acids) < CUTOFF for c in chgN):
                        new_hit += 1
                    cur = []
        rows.append({"receptor": rname, "ligand": lname, "poses": npose,
                     "old_any_N_pct": f"{100*old_hit/npose:.0f}",
                     "new_charged_N_pct": f"{100*new_hit/npose:.0f}"})

with open("docking/results_charged.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

LIG = ["L-Lys","L-Arg","L-Orn","D-Lys","L-Phe","L-Leu","L-Met","L-Trp"]
show = [r for r in ["inktoxin1","escapin","aplysianinA","snake_xtal"]
        if any(x["receptor"] == r for x in rows)]
print(f"\nsalt-bridge fraction (%), old metric -> new metric")
print(f"{'ligand':8s}" + "".join(f"{r[:13]:>16s}" for r in show))
for lg in LIG:
    line = f"{lg:8s}"
    for rc in show:
        m = [x for x in rows if x["receptor"] == rc and x["ligand"] == lg]
        line += f"{m[0]['old_any_N_pct']:>7s} ->{m[0]['new_charged_N_pct']:>5s}" if m else f"{'-':>16}"
    print(line)
print(f"\nwrote docking/results_charged.csv ({len(rows)} rows)")
