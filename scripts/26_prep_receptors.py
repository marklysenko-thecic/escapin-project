"""Phase 5.3a: prepare the four AlphaFold models as docking receptors.

Each AF3 model holds protein + FAD in one mmCIF. Meeko needs PDB input and
needs FAD kept (it shapes the pocket and its phosphates are charged), so the
cif is converted with FAD retained and any sugars dropped.

The snake crystal receptor already exists as docking/control/receptor_mkA.pdbqt.
"""
import glob, subprocess, shutil, sys
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO, Select

OUT = Path("docking/receptors"); OUT.mkdir(parents=True, exist_ok=True)
MODELS = {
    "inktoxin1":   "fold_inktoxin1_apunctata_fad",
    "escapin":     "fold_escapin_acalifornica_fad",
    "aplysianinA": "fold_aplysianina_acalifornica_fad",
    "snake_model": "fold_laao_crhodostoma_fad",
}
if shutil.which("mk_prepare_receptor.py") is None and shutil.which("mk_prepare_receptor") is None:
    sys.exit("meeko receptor tool not found - run: conda activate dock")
MK = "mk_prepare_receptor.py" if shutil.which("mk_prepare_receptor.py") else "mk_prepare_receptor"

class Keep(Select):
    def accept_residue(self, r):
        return r.id[0] == " " or r.get_resname() == "FAD"

cif, io = MMCIFParser(QUIET=True), PDBIO()
for name, folder in MODELS.items():
    hits = glob.glob(f"models/af3/{folder}/*model_0.cif")
    if not hits:
        print(f"{name}: model_0.cif not found"); continue
    st = cif.get_structure(name, hits[0])
    pdb = OUT / f"{name}.pdb"
    io.set_structure(st); io.save(str(pdb), Keep())
    nat = sum(1 for l in open(pdb) if l.startswith(("ATOM", "HETATM")))
    has_fad = any(l[17:20].strip() == "FAD" for l in open(pdb) if l.startswith("HETATM"))
    print(f"{name}: {nat} atoms written, FAD present: {has_fad}")

    r = subprocess.run([MK, "--read_pdb", str(pdb), "-o", str(OUT / name), "-p",
                        "--default_altloc", "A"], capture_output=True, text=True)
    pq = OUT / f"{name}.pdbqt"
    if not pq.exists():
        print(f"   FAILED - last lines:\n   " + "\n   ".join(
            (r.stdout + r.stderr).strip().splitlines()[-6:]))
        continue
    q = sum(float(l.split()[-2]) for l in open(pq) if l.startswith(("ATOM", "HETATM")))
    n = sum(1 for l in open(pq) if l.startswith(("ATOM", "HETATM")))
    print(f"   {pq.name}: {n} atoms, net charge {q:+.2f}")
