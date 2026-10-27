# Structural basis of substrate specificity in sea-hare L-amino acid oxidases

Escapin, from the ink of *Aplysia californica*, oxidises only L-lysine and
L-arginine. Related snake-venom enzymes prefer hydrophobic substrates. This
project asks why, using published data and public tools only.

## Findings

- Sea-hare pockets carry acidic residues (Glu207, Asp220, Glu224) where the
  snake enzyme has hydrophobic ones
- Pocket electrostatic potential: ink toxin 1 −18.6 kT/e, snake +8.3
  (crystal 2IID); five AlphaFold models per protein
- Those residues are conserved within *Aplysia* (9/10 sequences) but vary
  across other molluscs (13 different residues at position 207)
- One *A. californica* protein lacks them. Its pocket is the most positive
  measured (+13.6), and docking agrees, predicting it does not prefer
  lysine. Untested experimentally.
- A single Ser→Gly substitution beside the catalytic glutamate shifts
  escapin's predicted pocket potential by 13.3 kT/e

## Reproducing

Scripts are numbered in run order. Two conda environments: `escapin` for
sequence and structure analysis, `dock` for docking.

AlphaFold models are not included (size). Regenerate them by submitting the
sequences in `models/*_mature.fasta` to AlphaFold Server with FAD as a
ligand, one copy each, then unzip into `models/af3/`.

`notes/` holds the working log, including errors found and how.

## Data

| File | Contents |
|---|---|
| `results/apbs_ens/pocket_potential_ensemble.csv` | Pocket potentials, 5 models × 4 proteins |
| `results/apbs_outlier/outlier_potential.csv` | The outlier prediction test |
| `docking/results_poses.csv` | 1050 docked poses with metrics |
| `docking/outlier_salt.csv` | Docking cross-check of the prediction |
| `results/dataset_metadata.csv` | Sequence curation, 286 entries |
| `results/consurf/grades.txt` | Conservation grades |

## Tools

MAFFT, trimAl, IQ-TREE 3, CD-HIT, SignalP 6.0, AlphaFold Server, US-align,
PDB2PQR/PROPKA, APBS, Meeko, AutoDock Vina, RDKit, ConSurf, Biopython.
