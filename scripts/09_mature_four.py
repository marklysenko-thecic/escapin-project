"""Cut signal peptides from the four modelling targets."""
from Bio import SeqIO
from pathlib import Path

# REPLACE each 0 with the first number from SignalP's cleavage site
CLEAVAGE = {
    "escapin_Acalifornica":      18,
    "aplysianinA_Acalifornica":   19,
    "inktoxin1_Apunctata":        18,
    "LAAO_Crhodostoma":           19,
}

for name, pos in CLEAVAGE.items():
    f = Path(f"models/{name}_full.fasta")
    if not f.exists():
        print(f"MISSING {f}"); continue
    r = next(SeqIO.parse(f, "fasta"))
    mature = r.seq[pos:]
    out = Path(f"models/{name}_mature.fasta")
    out.write_text(f">{name}_mature signal peptide 1-{pos} removed\n{mature}\n")
    print(f"{name:28s} {len(r.seq)} -> {len(mature)} aa, starts {mature[:6]}")
