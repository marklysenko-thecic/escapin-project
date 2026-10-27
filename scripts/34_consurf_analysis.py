"""Phase 7.4: read the ConSurf grades and test H3.

H3 claims the substrate-pocket residues are conserved across sea-hare LAAOs
and their mollusc relatives. ConSurf grades each residue 1 (variable) to
9 (conserved) from evolutionary rates over the 37-sequence mollusc alignment.

The test is whether pocket residues score higher than the rest of the protein.
A Mann-Whitney U with a permutation p-value is used rather than a t-test,
because grades are ordinal and bounded.

ConSurf marks low-confidence grades with an asterisk; those are reported
separately, since a grade of 9 based on thin data is weaker evidence.
"""
import re, sys
from pathlib import Path
import numpy as np

SRC = Path("results/consurf/grades.txt")
POCKET = {75:"M",76:"R",103:"F",202:"H",203:"V",204:"F",205:"T",206:"S",
          207:"E",220:"D",224:"E",369:"Q",388:"Y",460:"W",499:"S",500:"W"}
ACIDIC = {207, 220, 224}

if not SRC.exists():
    sys.exit(f"missing {SRC}")

rows = {}
for line in open(SRC, errors="ignore"):
    # typical row: POS  SEQ  3LATOM  SCORE  COLOR  B/E  FUNCTION  MSA DATA ...
    m = re.match(r"\s*(\d+)\s+([A-Z])\s+(\S+)\s+(-?\d+\.\d+)\s+(\d)(\*?)", line)
    if m:
        pos, aa, score, grade, low = (int(m.group(1)), m.group(2),
                                      float(m.group(4)), int(m.group(5)), bool(m.group(6)))
        rows[pos] = {"aa": aa, "score": score, "grade": grade, "low_conf": low}

if not rows:
    print("could not parse; first 15 lines of the file:\n")
    for i, l in enumerate(open(SRC, errors="ignore")):
        if i >= 15: break
        print("   " + l.rstrip())
    sys.exit()

print(f"parsed {len(rows)} residues\n")

print(f"{'pos':>5} {'aa':>3} {'ConSurf aa':>11} {'score':>8} {'grade':>6} {'conf':>6}")
missing, mismatch = [], []
for p, expect in sorted(POCKET.items()):
    r = rows.get(p)
    if r is None:
        missing.append(p); continue
    if r["aa"] != expect: mismatch.append((p, expect, r["aa"]))
    tag = "  <-- acidic" if p in ACIDIC else ""
    print(f"{p:>5} {expect:>3} {r['aa']:>11} {r['score']:8.3f} {r['grade']:>6} "
          f"{'low' if r['low_conf'] else 'ok':>6}{tag}")

if missing:  print(f"\nnot found in the file: {missing}")
if mismatch:
    print(f"\nRESIDUE MISMATCH - numbering may be offset:")
    for p, e, g in mismatch: print(f"   position {p}: expected {e}, ConSurf has {g}")
    print("   (check whether a +18 precursor offset applies)")

pock = [rows[p]["grade"] for p in POCKET if p in rows]
rest = [v["grade"] for p, v in rows.items() if p not in POCKET]
print(f"\nmean grade: pocket {np.mean(pock):.2f} (n={len(pock)}), "
      f"rest of protein {np.mean(rest):.2f} (n={len(rest)})")
print(f"pocket at grade 9: {sum(1 for g in pock if g==9)}/{len(pock)}; "
      f"rest at grade 9: {100*sum(1 for g in rest if g==9)/len(rest):.0f}%")

# Mann-Whitney U as the probability a random pocket residue outranks a random other
u = sum((a > b) + 0.5*(a == b) for a in pock for b in rest) / (len(pock)*len(rest))
rng = np.random.default_rng(0)
allg = np.array(pock + rest)
null = []
for _ in range(20000):
    s = rng.permutation(allg)
    a, b = s[:len(pock)], s[len(pock):]
    null.append(sum((x > y) + 0.5*(x == y) for x in a for y in b) / (len(a)*len(b)))
p = (np.sum(np.array(null) >= u) + 1) / (len(null) + 1)
print(f"\nP(pocket residue more conserved than a random other) = {u:.3f}")
print(f"permutation p = {p:.4f}  (20000 shuffles, one-sided)")
print("  0.5 would mean no difference")
