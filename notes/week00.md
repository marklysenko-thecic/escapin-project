# Week 0 — Setup
Date:3.10.2026
## Done
-- Escapin from UniProt: Q6IWZ0, 535 aa, reviewed
- SignalP 6.0: cleavage between 18 and 19 Probability 0.979813 mature protein 517 aa
- NCBI search aplysianin[Title]: 3 results

## Findings
- aplysianin A (AAN78211.1, 556 aa): albumen gland (Cummins et al. 2004 JBC)
- Escapin: ink gland
- Paralogs: 61% identity, 100% cover → H4: different pockets/specificity?
- Escapin vs dactylomelin-P: alignment starts at escapin residue 19 (Query 19 = Sbjct 1), 25/26 identity → SignalP cleavage 18/19 confirmed by Edman data from a related species. Difference: R23 (escapin) / S (dactylomelin). Start of the FAD-binding motif visible (GxGxxG).
  
## Open questions
- Did anyone measure aplysianin A's substrates?
  
## Tool versions
(base) user@WSMSK1ML158 ~ % conda activate escapin
(escapin) user@WSMSK1ML158 ~ % mafft --version
trimal --version
iqtree2 --version || iqtree3 --version
cd-hit -h | head -3
blastp -version
python -c "import Bio; print(Bio.__version__)"
v7.526 (2024/Apr/26)
trimAl v1.5.rev1 build[2025-11-25]
zsh: command not found: iqtree2
IQ-TREE version 3.1.4 for MacOS ARM 64-bit built Sep 29 2026
Developed by Bui Quang Minh, Thomas Wong, Nhan Ly-Trong, Huaiyan Ren
Contributed by Lam-Tung Nguyen, Dominik Schrempf, Chris Bielow,
Olga Chernomor, Michael Woodhams, Diep Thi Hoang, Heiko Schmidt
		====== CD-HIT version 4.8.1 (built on Apr 24 2025) ======
Usage: cd-hit [Options] 
blastp: 2.17.0+
 Package: blast 2.17.0, build Aug 22 2026 11:50:41
1.88


## Next step
-
