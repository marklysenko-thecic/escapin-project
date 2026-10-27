"""Show the alignment window around each mapped pocket column."""
from Bio import AlignIO

SNAKE = "sp|P81382|OXLA_CALRH"
IDS = [SNAKE, "sp|Q6IWZ0|OXLA_APLCA", "AAR14185.1", "NP_001191524.1"]
LAB = ["snake ", "escap ", "inktx ", "aplyA "]
POCKET = [89, 90, 208, 223, 227, 322, 372, 374, 430, 464, 465]

aln = {r.id: str(r.seq) for r in AlignIO.read("data/processed/aligned.fasta","fasta")}
col, n = {}, 0
for i, c in enumerate(aln[SNAKE]):
    if c != "-":
        n += 1; col[n] = i

for pos in POCKET:
    c = col[pos + 18]
    lo, hi = max(0, c-6), c+7
    print(f"\n--- snake {aln[SNAKE][c]}{pos}  (alignment column {c}) ---")
    for i, lab in zip(IDS, LAB):
        w = aln[i][lo:hi]
        mark = " " * (c-lo) + "^"
        print(f"  {lab} {w}")
    print(f"         {mark}")
