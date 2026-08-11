#!/usr/bin/env python3
"""Reassess GNAZ(A)-RGSZ1(B) interface contacts across rank_1/2/3.
H-bonds are taken from ChimeraX (addh + hbonds, saved to hb_rankN.txt) - i.e. proper
donor/acceptor + distance + angle geometry. Salt bridge / hydrophobic / vdW are computed
here. A contact is 'common' if present (min heavy-atom <=4.0 A, or hydrophobic C-C <=4.5 A)
in ALL three models. pLDDT>=85 (both residues, min over the 3 models) filter reproduces the
'filtered' set. Writes common_contacts_filtered.csv (replacing it)."""
import os, math, re, csv
DIR=r"C:\Users\richi\Documents\RGSZ1 Human (262-378)-GNAZ-GTP Complex 4contact"
RANKS=["rank_1","rank_2","rank_3"]
CUT=4.0; SALT=4.0; HPHOB=4.5; HUMAN_OFFSET=261
POS={("ARG",a) for a in ("NH1","NH2","NE")}|{("LYS","NZ")}|{("HIS","ND1"),("HIS","NE2")}
NEG={("ASP",a) for a in ("OD1","OD2")}|{("GLU",a) for a in ("OE1","OE2")}
HYD={"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","CYS","TYR"}
PRIO={"Salt bridge":0,"H-bond":1,"Hydrophobic":2,"van der Waals":3}

def parse_pdb(path):
    A={};B={}  # resnum -> (resname, [(atom, x,y,z)]) ; plddt per residue
    plddt={"A":{}, "B":{}}
    for l in open(path):
        if l[:6] not in ("ATOM  ","HETATM"): continue
        ch=l[21]
        if ch not in ("A","B"): continue
        nm=l[12:16].strip()
        if nm[0]=="H" or (len(nm)>1 and nm[0].isdigit() and nm[1]=="H"): continue
        rn=l[17:20].strip(); ri=int(l[22:26])
        xyz=(float(l[30:38]),float(l[38:46]),float(l[46:54])); b=float(l[60:66])
        d=A if ch=="A" else B
        d.setdefault(ri,[rn,[]])[1].append((nm,xyz))
        plddt[ch].setdefault(ri,[]).append(b)
    pl={ch:{r:sum(v)/len(v) for r,v in plddt[ch].items()} for ch in ("A","B")}
    return A,B,pl

def parse_hb(path):
    """ChimeraX saveFile -> set of (a_resnum, b_resnum) residue pairs with an H-bond."""
    pairs=set()
    for l in open(path):
        m=re.findall(r"#\d+/([AB])\s+([A-Z0-9]{2,3})\s+(\d+)\s+\S+",l)
        if len(m)<2: continue
        # first two matches = donor, acceptor (third, if any, = hydrogen on donor residue)
        (c1,_,r1),(c2,_,r2)=m[0],m[1]
        if {c1,c2}=={"A","B"}:
            a=int(r1) if c1=="A" else int(r2)
            b=int(r1) if c1=="B" else int(r2)
            pairs.add((a,b))
    return pairs

def pair_metrics(A,B):
    """for every A-res x B-res: min heavy dist, salt flag, hydrophobic(C-C<=4.5) flag."""
    out={}
    for ar,(arn,aatoms) in A.items():
        for br,(brn,batoms) in B.items():
            mind=9e9; salt=False; hyd=False
            for an,ax in aatoms:
                for bn,bx in batoms:
                    d=math.dist(ax,bx)
                    if d>=HPHOB+0.01 and d>mind:
                        pass
                    if d<mind: mind=d
                    if d<=SALT and (((arn,an) in POS and (brn,bn) in NEG) or ((arn,an) in NEG and (brn,bn) in POS)):
                        salt=True
                    if d<=HPHOB and an[0]=="C" and bn[0]=="C" and arn in HYD and brn in HYD:
                        hyd=True
            out[(ar,arn,br,brn)]=(mind,salt,hyd)
    return out

data={}
for tag in RANKS:
    A,B,pl=parse_pdb(os.path.join(DIR,tag+".pdb"))
    hb=parse_hb(os.path.join(DIR,"hb_"+tag.replace("_","")+".txt"))
    pm=pair_metrics(A,B)
    data[tag]=dict(pm=pm,hb=hb,pl=pl)

def classify(tag,key):
    mind,salt,hyd=data[tag]["pm"][key]
    ar,arn,br,brn=key
    present = (mind<=CUT) or hyd
    if not present: return None,mind
    if salt: return "Salt bridge",mind
    if (ar,br) in data[tag]["hb"]: return "H-bond",mind
    if hyd: return "Hydrophobic",mind
    if mind<=CUT: return "van der Waals",mind
    return None,mind

# universe = any pair present in any rank
allkeys=set()
for tag in RANKS:
    for key in data[tag]["pm"]:
        t,_=classify(tag,key)
        if t: allkeys.add(key)

rows=[]
for key in allkeys:
    types=[];dists=[]
    for tag in RANKS:
        t,d=classify(tag,key); types.append(t); dists.append(d)
    if any(t is None for t in types): continue          # not present in all 3 -> not common
    ar,arn,br,brn=key
    consensus=sorted(types,key=lambda t:PRIO[t])[0]
    consistent = "yes" if len(set(types))==1 else "no"
    plA=min(data[tag]["pl"]["A"][ar] for tag in RANKS)
    plB=min(data[tag]["pl"]["B"][br] for tag in RANKS)
    rows.append(dict(key=key,ar=ar,arn=arn,br=br,brn=brn,types=types,consensus=consensus,
                     consistent=consistent,dists=dists,plA=plA,plB=plB,
                     mean=sum(dists)/3))

# pLDDT>=85 both chains (min over models) -> 'filtered'
filt=[r for r in rows if r["plA"]>=85 and r["plB"]>=85]
filt.sort(key=lambda r:(r["ar"],r["br"]))

# ---- write CSV ----
OUT=os.path.join(DIR,"common_contacts_filtered.csv")
with open(OUT,"w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["GNAZ_residue","RGSZ1_residue_human","RGSZ1_residue_file","contact_type",
                "consistent_type","mindist_rank1","mindist_rank2","mindist_rank3","mean_dist"])
    for r in filt:
        w.writerow([f"{r['arn']}{r['ar']}", f"{r['brn']}{r['br']+HUMAN_OFFSET}", f"{r['brn']}{r['br']}",
                    r["consensus"], r["consistent"],
                    f"{r['dists'][0]:.2f}",f"{r['dists'][1]:.2f}",f"{r['dists'][2]:.2f}",f"{r['mean']:.2f}"])
print("wrote",OUT,"with",len(filt),"common filtered contacts\n")

# ---- reconciliation vs the previous 27 rows (hard-coded from old CSV) ----
OLD={ # (GNAZ, file_resnum) -> old_type
 ("TYR75",101):"van der Waals",("ARG83",71):"H-bond",("ASP180",101):"H-bond",("ASP180",102):"H-bond",
 ("MET181",98):"Hydrophobic",("MET181",102):"van der Waals",("THR182",102):"van der Waals",
 ("THR183",27):"H-bond",("THR183",102):"H-bond",("THR183",103):"van der Waals",("THR183",106):"van der Waals",
 ("GLY184",22):"van der Waals",("ILE185",22):"H-bond",("ILE185",23):"Hydrophobic",("VAL186",106):"van der Waals",
 ("SER207",65):"van der Waals",("SER207",66):"van der Waals",("SER207",67):"H-bond",("GLU208",67):"H-bond",
 ("LYS211",23):"H-bond",("LYS211",26):"Salt bridge",("HIS214",23):"van der Waals",("TYR236",69):"H-bond",
 ("TYR236",70):"H-bond",("TYR236",71):"H-bond",("GLU237",70):"van der Waals",("ASN239",70):"H-bond"}
new={(f"{r['arn']}{r['ar']}",r["br"]):r["consensus"] for r in filt}
print("=== RECONCILIATION vs previous CSV (27 rows) ===")
print("CHANGED TYPE:")
for k,ot in sorted(OLD.items()):
    if k in new and new[k]!=ot:
        print(f"  {k[0]:8s}-file{k[1]:<4d} {ot:14s} -> {new[k]}")
print("DROPPED (no longer common in all 3 / not present):")
for k,ot in sorted(OLD.items()):
    if k not in new: print(f"  {k[0]:8s}-file{k[1]:<4d} (was {ot})")
print("NEW rows not in previous CSV:")
for k in sorted(new):
    if k not in OLD: print(f"  {k[0]:8s}-file{k[1]:<4d} = {new[k]}")
print(f"\nprev=27  new={len(filt)}")
