"""Phase 7.5: clade-level conservation of every pocket residue.

ConSurf grades conservation across all 37 molluscs, so a residue conserved
only within Aplysia scores as variable. That is what happened to Glu207
(grade 1, yet E in 9 of 10 sea hares).

This measures conservation separately inside and outside the sea-hare clade,
directly from the alignment, with no rate model - appropriate because the
10 Aplysia sequences are 93-99% identical and carry too few substitutions
for rate estimation.

Reported per position:
  Aplysia  - the majority residue and how many of 10 sequences carry it
  others   - the same for the 27 non-Aplysia molluscs, plus how many
             distinct residues appear there
  pattern  - CLADE-SPECIFIC when conserved in Aplysia but variable outside
"""
from collections import Counter
from Bio import AlignIO
import numpy as np

POCKET = {75:"M", 76:"R", 103:"F", 202:"H", 203:"V", 204:"F", 205:"T",
          206:"S", 207:"E", 220:"D", 224:"E", 369:"Q", 388:"Y",
          460:"W", 499:"S", 500:"W"}
ACIDIC = {207, 220, 224}
OFFSET = 18          # mature numbering -> precursor numbering, verified in Phase 4
GRADE = {75:9, 76:9, 103:2, 202:4, 203:4, 204:8, 205:3, 206:6,
         207:1, 220:3, 224:4, 369:8, 388:9, 460:6, 499:4, 500:9}

aln = AlignIO.read("data/processed/aligned.fasta", "fasta")
q = next(r for r in aln if "Q6IWZ0" in r.id)

def species(r):
    if "OS=" in r.description:
        return r.description.split("OS=")[1].split(" OX=")[0]
    return r.description.split("[")[-1].rstrip("]")

apl = [r for r in aln if species(r).startswith("Aplysia ")]
oth = [r for r in aln if not species(r).startswith("Aplysia ")
       and "CALRH" not in r.id and "HUMAN" not in r.id]
print(f"{len(apl)} Aplysia, {len(oth)} other molluscs (vertebrate outgroups excluded)\n")

# map mature position -> alignment column
col, n = {}, 0
for i, c in enumerate(str(q.seq)):
    if c != "-":
        n += 1
        col[n - OFFSET] = i          # n is precursor; subtract to get mature

print(f"{'pos':>4} {'exp':>4} {'ConSurf':>8} | {'Aplysia':>14} | {'other molluscs':>22} | pattern")
print("-" * 86)
rows = []
for p in sorted(POCKET):
    i = col.get(p)
    if i is None:
        print(f"{p:>4}  column not found"); continue
    ca = Counter(str(r.seq[i]) for r in apl)
    co = Counter(str(r.seq[i]) for r in oth)
    ta, na = ca.most_common(1)[0]
    to, no = co.most_common(1)[0]
    fa, fo = na / len(apl), no / len(oth)
    ndist = len([k for k in co if k != "-"])
    if fa >= 0.9 and fo < 0.5:
        pat = "CLADE-SPECIFIC"
    elif fa >= 0.9 and fo >= 0.9:
        pat = "conserved in both"
    elif fa < 0.9:
        pat = "variable in Aplysia too"
    else:
        pat = "partly shared"
    star = " *" if p in ACIDIC else "  "
    print(f"{p:>4} {POCKET[p]:>4} {GRADE.get(p,'?'):>8} | {ta} in {na:>2}/10 ({fa:4.0%}) | "
          f"{to} in {no:>2}/{len(oth)} ({fo:4.0%}), {ndist:>2} kinds | {pat}{star}")
    rows.append((p, fa, fo, pat))

print("\n* = acidic residue central to H1")
cs = [r for r in rows if r[3] == "CLADE-SPECIFIC"]
print(f"\nclade-specific positions: {len(cs)} of {len(rows)} -> "
      f"{', '.join(f'{POCKET[p]}{p}' for p, _, _, _ in cs) if cs else 'none'}")

# does ConSurf grade track conservation inside Aplysia, or outside it?
g = np.array([GRADE[p] for p, _, _, _ in rows], float)
fa = np.array([r[1] for r in rows]); fo = np.array([r[2] for r in rows])
def corr(a, b):
    a, b = a - a.mean(), b - b.mean()
    return float(a @ b / (np.sqrt((a**2).sum() * (b**2).sum()) + 1e-12))
print(f"\nConSurf grade vs conservation INSIDE Aplysia:  r = {corr(g, fa):+.2f}")
print(f"ConSurf grade vs conservation OUTSIDE Aplysia: r = {corr(g, fo):+.2f}")
print("a higher second value means the grade is driven by the non-Aplysia majority")
