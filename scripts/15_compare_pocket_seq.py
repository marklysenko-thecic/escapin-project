"""Phase 4.5: is escapin's pocket region identical to ink toxin 1's?

If the sequences match but the models differ in confidence, the low pLDDT
is modelling noise, not a real structural difference.
"""
from Bio import SeqIO

def load(p):
    return str(next(SeqIO.parse(p, "fasta")).seq)

esc = load("models/escapin_Acalifornica_mature.fasta")
ink = load("models/inktoxin1_Apunctata_mature.fasta")
apl = load("models/aplysianinA_Acalifornica_mature.fasta")

POS = [75, 76, 103, 202, 203, 204, 205, 206, 207, 220, 224, 341, 369, 388, 390, 460, 496, 499, 500]

print(f"escapin {len(esc)} aa, inktoxin1 {len(ink)} aa, aplysianinA {len(apl)} aa\n")

ident = sum(a == b for a, b in zip(esc, ink))
print(f"escapin vs inktoxin1, position by position: {ident}/{min(len(esc),len(ink))} "
      f"= {100*ident/min(len(esc),len(ink)):.1f}% identical\n")

print(f"{'pos':>5} {'escapin':>8} {'inktox1':>8} {'same?':>6}")
for p in POS:
    a = esc[p-1] if p <= len(esc) else "-"
    b = ink[p-1] if p <= len(ink) else "-"
    print(f"{p:>5} {a:>8} {b:>8} {'yes' if a == b else 'NO':>6}")

print("\nsequence window 198-230:")
print(f"  escapin    {esc[197:230]}")
print(f"  inktoxin1  {ink[197:230]}")
print(f"  position   198{' '*27}230")
