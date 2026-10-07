"""Step 1.4c: assemble the analysis set.

Rules (stated for the methods section):
  1. ALL full-length Aplysia sequences are kept, unclustered.
  2. Other molluscs: one representative per genus, chosen by a quality score.
  3. Sequences flagged partial / LOW QUALITY, or outside 430-650 aa, are excluded.
  4. Exact duplicates are removed.

Quality score for picking a genus representative (higher is better):
  +3  curated RefSeq (NP_) or Swiss-Prot
  +2  GenBank submission (direct author deposit)
  +1  length within 20 aa of escapin's 535
  -1  "-like" in the name (predicted by similarity only)
  -1  "isoform X" (one of several predicted splice variants)
"""
from Bio import SeqIO
from pathlib import Path

DROP = {"NP_001191570.1", "APIT_Aplysia_punctata"}   # verified duplicates
ESCAPIN_LEN = 535

def species(r):
    if "OS=" in r.description:
        return r.description.split("OS=")[1].split(" OX=")[0]
    return r.description.split("[")[-1].rstrip("]")

def usable(r):
    d = r.description.lower()
    if "partial" in d or "low quality" in d:
        return False
    return 430 <= len(r.seq) <= 650

def score(r):
    s, d = 0, r.description.lower()
    if r.id.startswith(("NP_", "sp|")) or "|" in r.id:
        s += 3
    elif not r.id.startswith("XP_"):
        s += 2
    if abs(len(r.seq) - ESCAPIN_LEN) <= 20:
        s += 1
    if "-like" in d:
        s -= 1
    if "isoform x" in d:
        s -= 1
    return s

seen, apl = set(), []
for p in ["data/raw/escapin_uniprot.fasta", "data/raw/blast_mollusca.fasta"]:
    for r in SeqIO.parse(p, "fasta"):
        if not species(r).startswith("Aplysia "):
            continue
        if r.id in DROP or str(r.seq) in seen or not usable(r):
            continue
        apl.append(r); seen.add(str(r.seq))

best = {}
for r in SeqIO.parse("data/processed/nr90.fasta", "fasta"):
    sp = species(r)
    if sp.startswith("Aplysia ") or str(r.seq) in seen or not usable(r):
        continue
    genus = sp.split()[0]
    if genus not in best or score(r) > score(best[genus]):
        best[genus] = r

oth = sorted(best.values(), key=lambda r: species(r))

Path("data/processed").mkdir(parents=True, exist_ok=True)
SeqIO.write(apl + oth, "data/processed/analysis_set.fasta", "fasta")

print(f"Aplysia {len(apl)} + other genera {len(oth)} = {len(apl)+len(oth)} sequences\n")
print("Representative per genus:")
for r in oth:
    print(f"  {r.id:22s} {len(r.seq):4d} aa  score {score(r):+d}  {species(r)}")
