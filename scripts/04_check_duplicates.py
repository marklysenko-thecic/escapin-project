"""Check which Aplysia sequences are identical to each other."""
from Bio import SeqIO
from collections import defaultdict

seqs = {}
for p in ["data/raw/escapin_uniprot.fasta",
          "data/raw/blast_mollusca.fasta",
          "data/raw/apunctata_laao_patent.fasta"]:
    for r in SeqIO.parse(p, "fasta"):
        if "Aplysia" in r.description or "APIT" in r.id:
            seqs.setdefault(r.id, str(r.seq))

groups = defaultdict(list)
for k, v in seqs.items():
    groups[v].append(k)

print("IDENTICAL GROUPS:")
for v, ks in groups.items():
    if len(ks) > 1:
        print(f"  {len(v)} aa: {', '.join(ks)}")

print(f"\n{len(seqs)} sequences, {len(groups)} unique")
