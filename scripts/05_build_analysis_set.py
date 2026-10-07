"""Analysis set: ALL full-length Aplysia (unclustered) + thinned other molluscs."""
from Bio import SeqIO

DROP = {"NP_001191570.1", "APIT_Aplysia_punctata"}

def ok(r):
    d = r.description.lower()
    return not ("partial" in d or "low quality" in d) and 430 <= len(r.seq) <= 650

def species(r):
    if "OS=" in r.description:
        return r.description.split("OS=")[1].split(" OX=")[0]
    return r.description.split("[")[-1].rstrip("]")

apl, seen = [], set()
for p in ["data/raw/escapin_uniprot.fasta", "data/raw/blast_mollusca.fasta"]:
    for r in SeqIO.parse(p, "fasta"):
        s = species(r)
        if not s.startswith("Aplysia "): continue
        if r.id in DROP or str(r.seq) in seen or not ok(r): continue
        apl.append(r); seen.add(str(r.seq))

oth = []
for r in SeqIO.parse("data/processed/nr90.fasta", "fasta"):
    if species(r).startswith("Aplysia") or str(r.seq) in seen or not ok(r): continue
    oth.append(r); seen.add(str(r.seq))

print(f"Aplysia: {len(apl)}, other molluscs: {len(oth)}\n")
for r in apl:
    print(f"  {r.id:24s} {len(r.seq)} aa  {species(r)}")
