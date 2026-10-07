"""Step 1.1: download escapin from UniProt."""
import requests
from pathlib import Path

OUT = Path("data/raw")
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://rest.uniprot.org/uniprotkb/search"
QUERY = '(protein_name:escapin) AND (organism_name:"Aplysia californica")'

# 1) Metadata table: check what was found
meta = requests.get(BASE, params={
    "query": QUERY,
    "format": "tsv",
    "fields": "accession,id,protein_name,organism_name,length,reviewed",
})
meta.raise_for_status()
print(meta.text)
(OUT / "escapin_uniprot.tsv").write_text(meta.text)

# 2) The sequence itself in FASTA
fasta = requests.get(BASE, params={"query": QUERY, "format": "fasta"})
fasta.raise_for_status()
(OUT / "escapin_uniprot.fasta").write_text(fasta.text)
print(fasta.text[:300])