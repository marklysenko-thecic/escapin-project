"""Summarise AlphaFold 3 confidence for each job (top-ranked model_0)."""
import json, glob
from pathlib import Path
from Bio.PDB import MMCIFParser

p = MMCIFParser(QUIET=True)
for d in sorted(Path("models/af3").iterdir()):
    if not d.is_dir():
        continue
    summ = glob.glob(str(d / "*summary_confidences_0.json"))
    cif  = glob.glob(str(d / "*model_0.cif"))
    if not summ or not cif:
        print(f"{d.name}: files missing"); continue
    s = json.load(open(summ[0]))
    st = p.get_structure("m", cif[0])
    pl = [a.get_bfactor() for a in st.get_atoms() if a.get_id() == "CA"]
    hi = sum(x > 70 for x in pl) / len(pl)
    print(f"{d.name}")
    print(f"   pTM {s.get('ptm')}   ipTM {s.get('iptm')}   clash {s.get('has_clash')}")
    print(f"   mean pLDDT {sum(pl)/len(pl):.1f}   residues >70: {100*hi:.0f}%   n={len(pl)}")
