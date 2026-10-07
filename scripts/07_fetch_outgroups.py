"""Step 1.5: fetch outgroup sequences for rooting the tree.

P81382 — snake LAAO (Calloselasma rhodostoma), the PDB 2IID protein.
         Also the structural control for Phase 3 and docking control for Phase 5.
Q96RQ9 — human IL4I1, a vertebrate LAAO; second outgroup to test the root.
"""
import requests
from pathlib import Path

ACC = ["P81382", "Q96RQ9"]
out = []
for a in ACC:
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{a}.fasta")
    r.raise_for_status()
    out.append(r.text.strip())
    header = r.text.split("\n")[0]
    seq_len = len("".join(r.text.split("\n")[1:]))
    print(f"{a}: {seq_len} aa")
    print(f"   {header[:100]}")

Path("data/raw").mkdir(parents=True, exist_ok=True)
Path("data/raw/outgroups.fasta").write_text("\n".join(out) + "\n")
