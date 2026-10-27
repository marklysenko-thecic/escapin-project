"""Phase 4.6a: convert AF3 models to PDB and run PDB2PQR at pH 5 and pH 8.

AF3 writes mmCIF with the protein and FAD in one entity set. PDB2PQR needs
PDB input and handles only standard residues, so FAD is stripped here and
added back for the pocket measurement in the next script.
"""
import glob, subprocess, sys
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO, Select

MODELS = {
    "escapin":     "models/af3/fold_escapin_acalifornica_fad/*model_0.cif",
    "inktoxin1":   "models/af3/fold_inktoxin1_apunctata_fad/*model_0.cif",
    "aplysianinA": "models/af3/fold_aplysianina_acalifornica_fad/*model_0.cif",
    "snake":       "models/af3/fold_laao_crhodostoma_fad/*model_0.cif",
}
PHS = [5.0, 8.0]
OUT = Path("results/apbs")
OUT.mkdir(parents=True, exist_ok=True)

class ProteinOnly(Select):
    def accept_residue(self, r):
        return r.id[0] == " "

cif = MMCIFParser(QUIET=True)
io = PDBIO()

for name, pat in MODELS.items():
    hits = glob.glob(pat)
    if not hits:
        print(f"{name}: model not found at {pat}"); continue
    st = cif.get_structure(name, hits[0])
    pdb = OUT / f"{name}.pdb"
    io.set_structure(st)
    io.save(str(pdb), ProteinOnly())
    print(f"{name}: wrote {pdb}")

    for ph in PHS:
        tag = f"{name}_pH{int(ph)}"
        cmd = ["pdb2pqr", "--ff=AMBER", "--titration-state-method=propka",
               f"--with-ph={ph}", "--apbs-input", str(OUT / f"{tag}.in"),
               str(pdb), str(OUT / f"{tag}.pqr")]
        r = subprocess.run(cmd, capture_output=True, text=True)
        ok = (OUT / f"{tag}.pqr").exists()
        print(f"   pH {ph}: {'ok' if ok else 'FAILED'}")
        if not ok:
            print("   " + r.stderr.strip().split("\n")[-1][:160])
