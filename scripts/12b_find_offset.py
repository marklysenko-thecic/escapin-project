"""Find the offset between crystal numbering and alignment numbering."""
from Bio import AlignIO

SNAKE = "sp|P81382|OXLA_CALRH"
POCKET = {89:"M", 90:"R", 208:"N", 223:"H", 227:"F", 322:"R",
          372:"Y", 374:"I", 430:"I", 464:"G", 465:"W"}

aln = {r.id: str(r.seq) for r in AlignIO.read("data/processed/aligned.fasta","fasta")}
snake = aln[SNAKE]
col, n = {}, 0
for i, c in enumerate(snake):
    if c != "-":
        n += 1
        col[n] = i

for off in range(14, 26):
    hits = sum(1 for p, aa in POCKET.items()
               if col.get(p+off) is not None and snake[col[p+off]] == aa)
    mark = "  <<<" if hits == len(POCKET) else ""
    print(f"offset +{off}: {hits}/{len(POCKET)} match{mark}")
