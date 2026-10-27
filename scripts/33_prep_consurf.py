"""Phase 7.3: prepare the ConSurf submission for H3.

ConSurf grades each residue 1 (variable) to 9 (conserved) from evolutionary
rates inferred over an alignment. H3 is a claim about sea hares and their
mollusc relatives, so the vertebrate outgroups are removed: at ~19% identity
they would dominate the rate estimates for positions that vary only within
the group of interest.

Outputs, all under results/consurf/:
  alignment_molluscs.fasta  the 37 mollusc sequences, gap-only columns removed
  escapin_query.fasta       escapin alone, as ConSurf's query sequence
  escapin_chainA.pdb        the model, single chain, protein only
  pocket_positions.txt      the residues to look up in the result
"""
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

OUTGROUPS = {"sp|P81382|OXLA_CALRH", "sp|Q96RQ9|OXLA_HUMAN"}
QUERY = "sp|Q6IWZ0|OXLA_APLCA"
POCKET = [75, 76, 103, 202, 203, 204, 205, 206, 207, 220, 224, 369, 388, 460, 499, 500]
OUT = Path("results/consurf"); OUT.mkdir(parents=True, exist_ok=True)

aln = AlignIO.read("data/processed/aligned.fasta", "fasta")
print(f"full alignment: {len(aln)} sequences, {aln.get_alignment_length()} columns")

kept = MultipleSeqAlignment([r for r in aln if r.id not in OUTGROUPS])
print(f"outgroups removed: {len(kept)} sequences remain")
if len(kept) != len(aln) - 2:
    print("  WARNING: expected to drop exactly 2 - check the outgroup ids")

# drop columns that are now gaps in every remaining sequence
keep_cols = [i for i in range(kept.get_alignment_length())
             if any(r.seq[i] != "-" for r in kept)]
trimmed = kept[:, keep_cols[0]:keep_cols[0]+1]
for i in keep_cols[1:]:
    trimmed += kept[:, i:i+1]
print(f"gap-only columns dropped: {kept.get_alignment_length()} -> {trimmed.get_alignment_length()}")

AlignIO.write(trimmed, OUT / "alignment_molluscs.fasta", "fasta")

q = next((r for r in trimmed if r.id == QUERY), None)
if q is None:
    raise SystemExit(f"query {QUERY} not found in the alignment")
seq = str(q.seq).replace("-", "")
(OUT / "escapin_query.fasta").write_text(f">escapin_Acalifornica\n{seq}\n")
print(f"query sequence written: {len(seq)} aa (expect 517)")

# single-chain protein-only PDB for ConSurf
src = Path("docking/receptors/escapin.pdb")
if src.exists():
    lines = [l for l in open(src) if l.startswith("ATOM")]
    (OUT / "escapin_chainA.pdb").write_text("".join(lines) + "END\n")
    print(f"structure written: {len(lines)} protein atoms (FAD excluded)")
else:
    print(f"NOTE: {src} not found - run scripts/26_prep_receptors.py first")

(OUT / "pocket_positions.txt").write_text(
    "Pocket residues to look up in the ConSurf output (mature numbering):\n" +
    "\n".join(f"  {p}" for p in POCKET) + "\n")
print(f"\npocket positions to check: {', '.join(map(str, POCKET))}")
print(f"all files in {OUT}/")
