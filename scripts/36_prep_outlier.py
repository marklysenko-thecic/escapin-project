"""Phase 7.6: prepare XP_005095149 and its Y207E mutant for AlphaFold.

XP_005095149 is the one Aplysia sequence with tyrosine, not glutamate, at the
key pocket position, and the one that branches off before both the ink-type
and albumen-type groups in the Phase 2 tree. It has no measured activity, so a
prediction about it is forward-looking rather than a retrodiction.

The natural sequence alone cannot isolate position 207 - it differs from
escapin at many positions. The Y207E mutant supplies the matched control, as
S206G did for escapin.

Position 207 is in MATURE numbering. This sequence has its own signal peptide,
so the cleavage site must be found for it rather than assumed from escapin.
"""
import re, sys
from pathlib import Path
from Bio import SeqIO, AlignIO

TARGET = "XP_005095149.2"
REF = "sp|Q6IWZ0|OXLA_APLCA"
ESC_MATURE_POS = 207          # escapin mature numbering
ESC_OFFSET = 18               # escapin precursor = mature + 18
OUT = Path("models"); OUT.mkdir(exist_ok=True)

# find the full-length sequence
rec = None
for p in ("data/raw/blast_mollusca.fasta", "data/raw/escapin_uniprot.fasta"):
    if not Path(p).exists():
        continue
    for r in SeqIO.parse(p, "fasta"):
        if r.id == TARGET:
            rec = r; break
    if rec: break
if rec is None:
    sys.exit(f"{TARGET} not found in data/raw/")
seq = str(rec.seq)
print(f"{TARGET}: {len(seq)} aa (precursor)")
print(f"  starts {seq[:25]}")

# locate the position equivalent to escapin mature 207, through the alignment
aln = {r.id: str(r.seq) for r in AlignIO.read("data/processed/aligned.fasta", "fasta")}
if TARGET not in aln or REF not in aln:
    sys.exit("target or reference missing from the alignment")
col, n = None, 0
for i, c in enumerate(aln[REF]):
    if c != "-":
        n += 1
        if n == ESC_MATURE_POS + ESC_OFFSET:
            col = i; break
tgt_pos, n = None, 0
for i, c in enumerate(aln[TARGET]):
    if c != "-":
        n += 1
    if i == col:
        tgt_pos = n; break
print(f"\nalignment column {col}: escapin '{aln[REF][col]}', {TARGET} '{aln[TARGET][col]}'")
print(f"  equivalent precursor position in {TARGET}: {tgt_pos}")
if aln[TARGET][col] != "Y":
    print(f"  WARNING: expected Y, found '{aln[TARGET][col]}'")
if seq[tgt_pos-1] != aln[TARGET][col]:
    sys.exit(f"  consistency check failed: seq has '{seq[tgt_pos-1]}' at {tgt_pos}")
print(f"  context: ...{seq[tgt_pos-6:tgt_pos-1]}[{seq[tgt_pos-1]}]{seq[tgt_pos:tgt_pos+5]}...")

(OUT / f"{TARGET}_full.fasta").write_text(f">{TARGET}\n{seq}\n")
print(f"\nwritten: models/{TARGET}_full.fasta")
print("\nNEXT: run this sequence through SignalP 6.0 (Eukarya) to get its own")
print("cleavage site, then rerun this script with CLEAVAGE set below.")

CLEAVAGE = 22      # <-- set to the first number of SignalP's 'CS pos: N-N+1'
if CLEAVAGE:
    mature = seq[CLEAVAGE:]
    mpos = tgt_pos - CLEAVAGE
    if mature[mpos-1] != "Y":
        sys.exit(f"after cleavage, position {mpos} is '{mature[mpos-1]}', expected Y")
    mut = mature[:mpos-1] + "E" + mature[mpos:]
    assert sum(a != b for a, b in zip(mature, mut)) == 1
    (OUT / "XP005095149_wt.fasta").write_text(f">XP005095149_wt\n{mature}\n")
    (OUT / "XP005095149_Y207E.fasta").write_text(f">XP005095149_Y207E\n{mut}\n")
    print(f"\nmature {len(mature)} aa, mutation at mature position {mpos}")
    print(f"  wt  ...{mature[mpos-6:mpos+5]}...")
    print(f"  mut ...{mut[mpos-6:mpos+5]}...")
    print("both written to models/ - submit each with FAD, Copies 1")
