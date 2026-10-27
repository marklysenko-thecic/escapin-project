"""Phase 4.9: build the position-206 swap mutants for AlphaFold.

escapin   S206G  -> should become ink-toxin-like if 206 controls the pocket
inktoxin1 G206S  -> should become escapin-like
Numbering is the mature protein (signal peptide removed), as in all pocket work.
"""
from Bio import SeqIO
from pathlib import Path

JOBS = [("escapin_Acalifornica", "S", "G", "escapin_S206G"),
        ("inktoxin1_Apunctata",  "G", "S", "inktoxin1_G206S")]
POS = 206

for src, wt, mut, out in JOBS:
    seq = str(next(SeqIO.parse(f"models/{src}_mature.fasta", "fasta")).seq)
    found = seq[POS - 1]
    if found != wt:
        raise SystemExit(f"{src}: expected {wt} at {POS}, found {found} - stopping")
    new = seq[:POS - 1] + mut + seq[POS:]
    assert len(new) == len(seq) and sum(a != b for a, b in zip(seq, new)) == 1
    Path(f"models/{out}.fasta").write_text(f">{out}\n{new}\n")
    print(f"{out}: {len(new)} aa, exactly 1 change")
    print(f"   wild type  ...{seq[POS-6:POS+5]}...")
    print(f"   mutant     ...{new[POS-6:POS+5]}...")
    print(f"               {' ' * 5}^ position {POS}\n")
