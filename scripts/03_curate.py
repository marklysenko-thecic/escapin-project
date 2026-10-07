"""Step 1.4: merge all sources, filter, build metadata table."""
import re, pandas as pd
from pathlib import Path
from Bio import SeqIO

SOURCES = [
    ("data/raw/escapin_uniprot.fasta", "UniProt"),
    ("data/raw/blast_mollusca.fasta", "BLAST_mollusca"),
    ("data/raw/apunctata_laao_patent.fasta", "patent_WO2004065415A2"),
]

rows, records = [], {}
for path, source in SOURCES:
    if not Path(path).exists():
        print(f"MISSING: {path}"); continue
    for r in SeqIO.parse(path, "fasta"):
        if r.id in records: continue
        records[r.id] = r
        desc = r.description
        sp = re.search(r"\[([^\]]+)\]", desc)
        species = sp.group(1) if sp else ("Aplysia californica" if source == "UniProt" else "Aplysia punctata")
        rows.append({
            "id": r.id,
            "species": species,
            "length": len(r.seq),
            "source": source,
            "partial": "partial" in desc.lower(),
            "low_quality": "low quality" in desc.lower(),
            "description": desc[:90],
        })

df = pd.DataFrame(rows)
df["keep"] = (~df.partial) & (~df.low_quality) & df.length.between(430, 650)
Path("results").mkdir(exist_ok=True)
df.to_csv("results/dataset_metadata.csv", index=False)

print(f"total {len(df)}, passing filters {df.keep.sum()}")
print(df[df.species.str.startswith("Aplysia")][["id","species","length","partial"]].to_string(index=False))
