"""Phase 3 prep: extract the four sequences to model."""
from Bio import SeqIO
from pathlib import Path

TARGETS = {
    "sp|Q6IWZ0|OXLA_APLCA": "escapin_Acalifornica",
    "NP_001191524.1":       "aplysianinA_Acalifornica",
    "AAR14185.1":           "inktoxin1_Apunctata",
    "sp|P81382|OXLA_CALRH": "LAAO_Crhodostoma",
}

found = {}
for p in ["data/raw/escapin_uniprot.fasta",
          "data/raw/blast_mollusca.fasta",
          "data/raw/outgroups.fasta"]:
    for r in SeqIO.parse(p, "fasta"):
        if r.id in TARGETS:
            found[r.id] = r

Path("models").mkdir(exist_ok=True)
for acc, name in TARGETS.items():
    if acc not in found:
        print(f"MISSING: {acc}"); continue
    r = found[acc]
    print(f"\n>{name}  ({acc}, {len(r.seq)} aa)")
    print(str(r.seq)[:60] + "...")
    (Path("models") / f"{name}_full.fasta").write_text(f">{name}\n{r.seq}\n")
