"""Phase 5.5: re-score docked poses using VERIFIED alpha-atom indices.

Reads docking/ligand_atom_index.json, written by the coordinate-matching step,
so the alpha carbon and alpha nitrogen are the real ones for each ligand.
The previous attempt assumed Meeko preserved RDKit's atom order; it does not,
and every geometric measurement was taken from the wrong atoms.

Reports each criterion separately before combining them, so a failing filter
is visible instead of silently zeroing the result.
"""
import csv, glob, json, re, sys
from pathlib import Path
import numpy as np

CA_N5_MAX, N_BASE_MAX, ANGLE_MAX, SALT_MAX = 5.0, 7.0, 75.0, 4.0
AD = {"A":"C","C":"C","OA":"O","O":"O","N":"N","NA":"N","SA":"S","S":"S"}
REC_PATH = {"inktoxin1":"docking/receptors/inktoxin1.pdbqt",
            "escapin":"docking/receptors/escapin.pdbqt",
            "aplysianinA":"docking/receptors/aplysianinA.pdbqt",
            "snake_model":"docking/receptors/snake_model.pdbqt",
            "snake_xtal":"docking/control/receptor_mkA.pdbqt"}
IDX = json.load(open("docking/ligand_atom_index.json"))

def rec_atoms(p):
    out=[]
    for l in open(p):
        if l.startswith(("ATOM","HETATM")):
            t=l.split()[-1].upper()
            out.append((l[17:20].strip(), l[12:16].strip(), AD.get(t,t[0]),
                        np.array([float(l[30+8*i:38+8*i]) for i in range(3)])))
    return out

cache={}
def info(name):
    if name not in cache:
        ra=rec_atoms(REC_PATH[name])
        n5=next((c for rn,an,el,c in ra if rn=="FAD" and an=="N5"), None)
        base=[c for rn,an,el,c in ra if rn=="HIS" and an in ("ND1","NE2")
              and n5 is not None and np.linalg.norm(c-n5)<12.0]
        acid=[c for rn,an,el,c in ra if rn in ("ASP","GLU")
              and an in ("OD1","OD2","OE1","OE2")
              and n5 is not None and np.linalg.norm(c-n5)<14.0]
        cache[name]=(n5,base,acid)
    return cache[name]

rows=[]
for f in sorted(glob.glob("docking/runs/*.pdbqt")):
    m=re.match(r"(.+?)__(.+?)__s(\d+)\.pdbqt", Path(f).name)
    if not m: continue
    rn,lg,seed=m.group(1),m.group(2),int(m.group(3))
    if rn not in REC_PATH or lg not in IDX: continue
    n5,bases,acids=info(rn)
    if n5 is None: continue
    iN,iCA,iC = IDX[lg]["N"], IDX[lg]["CA"], IDX[lg]["C"]

    poses,cur,scores=[],[],[]
    for l in open(f):
        if l.startswith("REMARK VINA RESULT"): scores.append(float(l.split()[3]))
        elif l.startswith(("ATOM","HETATM")): cur.append(l)
        elif l.startswith("ENDMDL"): poses.append(cur); cur=[]
    if cur: poses.append(cur)

    for i,(pl,sc) in enumerate(zip(poses,scores),1):
        h=[]
        for l in pl:
            t=l.split()[-1].upper()
            if t in ("HD","H"): continue
            h.append((AD.get(t,t[0]), np.array([float(l[30+8*i2:38+8*i2]) for i2 in range(3)])))
        if max(iN,iCA,iC)>=len(h): continue
        aN,aCA,aC = h[iN][1], h[iCA][1], h[iC][1]
        d_ca=float(np.linalg.norm(aCA-n5))
        d_base=min((float(np.linalg.norm(aN-b)) for b in bases), default=np.inf)
        # hydride direction: CA's fourth tetrahedral vertex, opposite its 3 heavy neighbours
        nb=np.array([c for j,(e,c) in enumerate(h) if j!=iCA])
        near=nb[np.argsort(np.linalg.norm(nb-aCA,axis=1))[:3]]
        hv=aCA-near.mean(axis=0); tv=n5-aCA
        ang=float(np.degrees(np.arccos(np.clip(hv@tv/(np.linalg.norm(hv)*np.linalg.norm(tv)+1e-9),-1,1))))
        salt=sum(1 for e,c in h if e=="N" and min((np.linalg.norm(c-a) for a in acids),default=99)<SALT_MAX)
        rows.append({"receptor":rn,"ligand":lg,"seed":seed,"pose":i,"score":f"{sc:.2f}",
                     "CA_to_N5":f"{d_ca:.2f}","N_to_base":f"{d_base:.2f}",
                     "angle":f"{ang:.0f}","salt":salt,
                     "c1_dist":int(d_ca<CA_N5_MAX),"c2_base":int(d_base<N_BASE_MAX),
                     "c3_angle":int(ang<ANGLE_MAX),
                     "catalytic":int(d_ca<CA_N5_MAX and d_base<N_BASE_MAX and ang<ANGLE_MAX)})

with open("docking/results_catalytic2.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"re-scored {len(rows)} poses\n")

LIG=["L-Lys","L-Arg","L-Orn","D-Lys","L-Phe","L-Leu","L-Met","L-Trp"]
RECS=["inktoxin1","escapin","aplysianinA","snake_model","snake_xtal"]

print("EACH CRITERION ALONE, pooled over receptors (% of poses passing)")
print(f"{'ligand':8s} {'CA<5A':>8} {'N-base<7A':>11} {'angle<75':>10} {'all three':>11} {'salt>0':>8}")
for lg in LIG:
    s=[r for r in rows if r["ligand"]==lg]
    print(f"{lg:8s} {100*sum(r['c1_dist'] for r in s)/len(s):8.0f}"
          f"{100*sum(r['c2_base'] for r in s)/len(s):11.0f}"
          f"{100*sum(r['c3_angle'] for r in s)/len(s):10.0f}"
          f"{100*sum(r['catalytic'] for r in s)/len(s):11.0f}"
          f"{100*sum(1 for r in s if r['salt']>0)/len(s):8.0f}")

print("\nCATALYTIC POSE FRACTION (%) by receptor")
print(f"{'ligand':8s}"+"".join(f"{r[:11]:>13s}" for r in RECS))
for lg in LIG:
    line=f"{lg:8s}"
    for rc in RECS:
        s=[r for r in rows if r["ligand"]==lg and r["receptor"]==rc]
        line += f"{100*sum(r['catalytic'] for r in s)/len(s):13.0f}" if s else f"{'-':>13}"
    print(line)

print("\nL vs D lysine (the stereochemistry control):")
for rc in RECS:
    a=[r for r in rows if r["receptor"]==rc and r["ligand"]=="L-Lys"]
    b=[r for r in rows if r["receptor"]==rc and r["ligand"]=="D-Lys"]
    if a and b:
        fa=100*sum(r["catalytic"] for r in a)/len(a); fb=100*sum(r["catalytic"] for r in b)/len(b)
        aa=np.mean([float(r["angle"]) for r in a]); ab=np.mean([float(r["angle"]) for r in b])
        print(f"  {rc:12s} L {fa:5.0f}%  D {fb:5.0f}%   mean angle L {aa:3.0f} deg  D {ab:3.0f} deg")
