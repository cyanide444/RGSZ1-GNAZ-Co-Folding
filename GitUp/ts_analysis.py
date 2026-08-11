#!/usr/bin/env python3
"""Transition-state comparison of the superimposed 1AGR (GNAI1.GDP.AlF4.Mg, TS mimic)
and rank_1 (GNAZ.GTP.Mg, ground state). Both PDBs are in the same (Matchmaker) frame."""
import math
P=r"C:\Users\richi\Documents\RGSZ1-GNAZ Project"
def parse(path):
    at=[]
    for l in open(path):
        if l[:6] not in ("ATOM  ","HETATM"): continue
        at.append(dict(rec=l[:6].strip(),name=l[12:16].strip(),resn=l[17:21].strip(),
            ch=l[21],ri=l[22:26].strip(),el=(l[76:78].strip() or l[12:16].strip()[0]),
            xyz=(float(l[30:38]),float(l[38:46]),float(l[46:54]))))
    return at
AG=parse(P+r"\_1agr_super.pdb")   # 1AGR (GNAI1=A, RGS4=E)
RK=parse(P+r"\_rank1_super.pdb")  # rank_1 (GNAZ=A, RGSZ1=B)
def d(a,b): return math.dist(a,b)
def get(at,ch=None,ri=None,resn=None,name=None):
    r=[a for a in at if (ch is None or a["ch"]==ch) and (ri is None or a["ri"]==str(ri))
       and (resn is None or a["resn"]==resn) and (name is None or a["name"]==name)]
    return r
def atom(at,**k):
    r=get(at,**k); return r[0]["xyz"] if r else None

# --- 1AGR TS-mimic anchors ---
AL=atom(AG,resn="ALF",name="AL")
GLN204={n:atom(AG,ch="A",ri=204,name=n) for n in ("CA","CB","CG","CD","OE1","NE2")}
# arginine finger = Arg whose guanidinium N is closest to AL
def arg_finger(at,target):
    best=None
    for a in at:
        if a["resn"]=="ARG" and a["name"] in ("NH1","NH2","NE"):
            dist=d(a["xyz"],target)
            if best is None or dist<best[0]: best=(dist,a["ch"],a["ri"],a["name"],a["xyz"])
    return best
AGfing=arg_finger([a for a in AG if a["ch"]=="A"],AL)
# catalytic water: HOH O closest to (AL and Gln204.NE2)
def cat_water(at,al,gln_ne2):
    best=None
    for a in at:
        if a["resn"]=="HOH":
            s=d(a["xyz"],al)+d(a["xyz"],gln_ne2)
            if best is None or s<best[0]: best=(s,a["xyz"],d(a["xyz"],al),d(a["xyz"],gln_ne2))
    return best
W=cat_water(AG,AL,GLN204["NE2"])
MGa=atom(AG,resn="MG",name="MG") or atom(AG,ch="A",name="MG")

# --- rank_1 ground-state anchors ---
GLN205={n:atom(RK,ch="A",ri=205,name=n) for n in ("CA","CB","CG","CD","OE1","NE2")}
# GTP gamma phosphate: among P atoms, classify by P-O-P bridges and ribose ester
Ps=[a for a in RK if a["resn"]=="LIG1" and a["el"]=="P"]
Os=[a for a in RK if a["resn"]=="LIG1" and a["el"]=="O"]
Cs=[a for a in RK if a["resn"]=="LIG1" and a["el"]=="C"]
def bonded(a,b,cut=1.9): return d(a["xyz"],b["xyz"])<cut
Pinfo={}
for p in Ps:
    boxy=[o for o in Os if bonded(p,o)]
    nbridge=sum(1 for o in boxy if sum(1 for q in Ps if q is not p and bonded(o,q))>=1)
    ester=any(sum(1 for c in Cs if bonded(o,c))>=1 for o in boxy if any(bonded(o,q) for q in Ps if q is not p))
    Pinfo[p["name"]]=dict(nbridge=nbridge,ester=ester,xyz=p["xyz"],boxy=[o["name"] for o in boxy])
# gamma = P with exactly 1 bridging O and NOT ester-linked to ribose
gamma=None;beta=None;alpha=None
for nm,i in Pinfo.items():
    if i["nbridge"]>=2: beta=nm
    elif i["ester"]: alpha=nm
    else: gamma=nm
PG=Pinfo[gamma]["xyz"] if gamma else None
# terminal O of gamma phosphate
GAMMA_O=[o["xyz"] for o in Os if bonded(o,Pinfo[gamma]) and o["name"]] if gamma else []
MGr=atom(RK,resn="LIG2") or [a["xyz"] for a in RK if a["el"]=="Mg"][0]
if MGr is None:
    MGr=[a["xyz"] for a in RK if a["el"]=="Mg"][0]
RKfing=arg_finger([a for a in RK if a["ch"]=="A"],PG)

print("="*70)
print("1AGR (GNAI1.GDP.AlF4.Mg = transition-state mimic)")
print("  AlF4 Al present:",AL is not None,"| Mg:",MGa is not None)
print("  arginine finger:", "Arg%s/%s @ %.2f A from Al"%(AGfing[2],AGfing[3],AGfing[0]))
print("  catalytic water: d(W-Al)=%.2f  d(W-Gln204.NE2)=%.2f"%(W[2],W[3]))
print("  Gln204.NE2 -> Al = %.2f   Gln204.OE1 -> Al = %.2f"%(d(GLN204["NE2"],AL),d(GLN204["OE1"],AL)))
print()
print("rank_1 (GNAZ.GTP.Mg = ground state)")
print("  GTP phosphates: alpha=%s beta=%s gamma=%s"%(alpha,beta,gamma))
print("  arginine finger:", "Arg%s/%s @ %.2f A from Pgamma"%(RKfing[2],RKfing[3],RKfing[0]))
print("  Gln205.NE2 -> Pgamma = %.2f   Gln205.OE1 -> Pgamma = %.2f"%(d(GLN205["NE2"],PG),d(GLN205["OE1"],PG)))
print()
print("="*70)
print("SUPERPOSITION-BASED COMPARISON (same frame)")
print("  Catalytic Gln overlay  Gln204(1AGR) vs Gln205(rank1):")
for n in ("CA","CB","CG","CD","OE1","NE2"):
    print("     %-4s  %.2f A"%(n,d(GLN204[n],GLN205[n])))
print("  gamma-phosphate mimic: AlF4-Al (1AGR) vs GTP-Pgamma (rank1) = %.2f A"%d(AL,PG))
print("  Mg(1AGR) vs Mg(rank1) = %.2f A"%d(MGa,MGr))
print("  arg-finger guanidinium: 1AGR-%s vs rank1-%s = %.2f A"%(AGfing[3],RKfing[3],d(AGfing[4],RKfing[4])))
print()
print("  TRANSITION-STATE PROBE - put 1AGR catalytic water into rank_1 frame:")
print("     d(W_1AGR -> rank1 Gln205.NE2) = %.2f  (vs %.2f in 1AGR to its own Gln)"%(d(W[1],GLN205["NE2"]),W[3]))
print("     d(W_1AGR -> rank1 GTP-Pgamma) = %.2f  (vs %.2f in 1AGR W->Al)"%(d(W[1],PG),W[2]))
print()
print("  local backbone fit (Ca of equivalent catalytic residues):")
for lbl,a4,r5 in [("catGln",204,205),("argFinger",int(AGfing[2]),int(RKfing[2])),
                  ("SwI-Thr",181,182)]:
    ca4=atom(AG,ch="A",ri=a4,name="CA"); ca5=atom(RK,ch="A",ri=r5,name="CA")
    if ca4 and ca5: print("     %-10s GNAI1 %d / GNAZ %d  Ca-Ca = %.2f A"%(lbl,a4,r5,d(ca4,ca5)))
