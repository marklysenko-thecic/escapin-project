"""Phase 4.2: map snake pocket positions onto the other three proteins."""
from Bio import AlignIO

SNAKE = "sp|P81382|OXLA_CALRH"
TARGETS = {"sp|Q6IWZ0|OXLA_APLCA": "escapin",
           "AAR14185.1": "inktoxin1", "NP_001191524.1": "aplysianinA"}
# crystal numbering (mature protein, starts at 1)
POCKET = [89, 90, 208, 223, 227, 322, 372, 374, 430, 464, 465]

aln = {r.id: str(r.seq) for r in AlignIO.read("data/processed/aligned.fasta", "fasta")}
snake = aln[SNAKE]

# map ungapped position -> alignment column
col = {}
n = 0
for i, c in enumerate(snake):
    if c != "-":
        n += 1
        col[n] = i

print(f"{'snake':>10} | " + " | ".join(f"{v:>11}" for v in TARGETS.values()))
print("-" * 54)
for pos in POCKET:
    # the crystal is the mature protein; the alignment uses the precursor (+19)
    c = col.get(pos + 18)
    if c is None:
        print(f"{pos:>10} | (not in alignment)"); continue
    row = f"{snake[c]}{pos:<9}"
    cells = [f"{aln[t][c]:>11}" for t in TARGETS]
    print(f"{row:>10} | " + " | ".join(cells))
