#!/usr/bin/env python3
"""Regenerate the GNAZ-RGSZ1 contact table directly from rank_1/2/3.cif (auth_seq_id
numbering, RGSZ1 chain B = 1-117). H-bonds from ChimeraX (hb_cif*.txt). Adds region
labels: GNAZ switch I/II/III/other; RGSZ1 Tesmer et al. (1997) helix/loop via the
RGS4 homology (RGS4 residue = RGSZ1 AF3 file + 61)."""
import os, math, re, csv
D=r"C:\Users\richi\Documents\RGSZ1 Human (262-378)-GNAZ-GTP Complex 4contact"
RANKS=["rank_1","rank_2","rank_3"]
CUT=4.0; SALT=4.0; HPHOB=4.5; HUMAN=261   # human RGS20 = AF3 file + 261
POS={("ARG",a) for a in ("NH1","NH2","NE")}|{("LYS","NZ")}|{("HIS","ND1"),("HIS","NE2")}
NEG={("ASP",a) for a in ("OD1","OD2")}|{("GLU",a) for a in ("OE1","OE2")}
HYD={"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","CYS","TYR"}
PRIO={"Salt bridge":0,"H-bond":1,"Hydrophobic":2,"van der Waals":3}

def parse_cif(path):
    A={};B={};pl={"A":{}, "B":{}}
    for l in open(path):
        if not (l.startswith("ATOM ") or l.startswith("HETATM ")): continue
        f=l.split()
        if len(f)<18: continue
        elem=f[2]; nm=f[3]; resn=f[5]; ch=f[15]
        if elem=="H": continue
        if ch not in ("A","B"): continue
        try:
            ri=int(f[7]); xyz=(float(f[10]),float(f[11]),float(f[12])); b=float(f[17])
        except ValueError: continue
        d=A if ch=="A" else B
        d.setdefault(ri,[resn,[]])[1].append((nm,xyz))
        pl[ch].setdefault(ri,[]).append(b)
    plm={ch:{r:sum(v)/len(v) for r,v in pl[ch].items()} for ch in ("A","B")}
    return A,B,plm

def parse_hb(path):
    pairs=set()
    for l in open(path):
        m=re.findall(r"#\d+/([AB])\s+([A-Z0-9]{2,3})\s+(\d+)\s+\S+",l)
        if len(m)<2: continue
        (c1,_,r1),(c2,_,r2)=m[0],m[1]
        if {c1,c2}=={"A","B"}:
            a=int(r1) if c1=="A" else int(r2); b=int(r1) if c1=="B" else int(r2)
            pairs.add((a,b))
    return pairs

def pair_metrics(A,B):
    out={}
    for ar,(arn,aa) in A.items():
        for br,(brn,ba) in B.items():
            mind=9e9; salt=False; hyd=False
            for an,ax in aa:
                for bn,bx in ba:
                    dd=math.dist(ax,bx)
                    if dd<mind: mind=dd
                    if dd<=SALT and (((arn,an) in POS and (brn,bn) in NEG) or ((arn,an) in NEG and (brn,bn) in POS)): salt=True
                    if dd<=HPHOB and an[0]=="C" and bn[0]=="C" and arn in HYD and brn in HYD: hyd=True
            out[(ar,arn,br,brn)]=(mind,salt,hyd)
    return out

# ---- region labels ----
def gnaz_region(n):
    if 175<=n<=188: return "Switch I"
    if 200<=n<=220: return "Switch II"
    if 228<=n<=239: return "Switch III"
    return "other"
def rgs4_element(n):   # Tesmer 1997 RGS4 boundaries (from 1AGR HELIX records)
    for a,b,name in [(53,61,"a1"),(62,62,"a1-a2 loop"),(63,68,"a2"),(69,69,"a2-a3 loop"),
        (70,82,"a3"),(83,85,"a3-a4 loop"),(86,100,"a4"),(101,105,"a4-a5 loop"),
        (106,118,"a5"),(119,130,"a5-a6 loop"),(131,142,"a6"),(143,149,"a6-a7 loop"),
        (150,162,"a7"),(163,163,"a7-a8 loop"),(164,170,"a8"),(171,171,"a8-a9 loop"),
        (172,175,"a9")]:
        if a<=n<=b: return name
    return "terminal/other"
def rgsz1_element(af3_file): return rgs4_element(af3_file+61)

data={}
for tag in RANKS:
    A,B,pl=parse_cif(os.path.join(D,tag+".cif"))
    hb=parse_hb(os.path.join(D,"hb_cif"+tag.split("_")[1]+".txt"))
    data[tag]=dict(pm=pair_metrics(A,B),hb=hb,pl=pl)

def classify(tag,key):
    mind,salt,hyd=data[tag]["pm"][key]; ar,arn,br,brn=key
    present=(mind<=CUT) or hyd
    if not present: return None,mind
    if salt: return "Salt bridge",mind
    if (ar,br) in data[tag]["hb"]: return "H-bond",mind
    if hyd: return "Hydrophobic",mind
    if mind<=CUT: return "van der Waals",mind
    return None,mind

allkeys=set()
for tag in RANKS:
    for key in data[tag]["pm"]:
        t,_=classify(tag,key)
        if t: allkeys.add(key)
rows=[]
for key in allkeys:
    types=[];dists=[]
    for tag in RANKS:
        t,dd=classify(tag,key); types.append(t); dists.append(dd)
    if any(t is None for t in types): continue
    ar,arn,br,brn=key
    plA=min(data[tag]["pl"]["A"][ar] for tag in RANKS)
    plB=min(data[tag]["pl"]["B"][br] for tag in RANKS)
    rows.append(dict(ar=ar,arn=arn,br=br,brn=brn,
        consensus=sorted(types,key=lambda t:PRIO[t])[0],
        consistent=("yes" if len(set(types))==1 else "no"),
        dists=dists,plA=plA,plB=plB,mean=sum(dists)/3))
filt=[r for r in rows if r["plA"]>=85 and r["plB"]>=85]
filt.sort(key=lambda r:(r["ar"],r["br"]))

OUT=os.path.join(D,"common_contacts_filtered.csv")
try:
    open(OUT,"a").close()
except PermissionError:
    OUT=os.path.join(D,"common_contacts_filtered_regen.csv")
    print("** original CSV locked (open in Excel?); writing to",os.path.basename(OUT))
with open(OUT,"w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f)
    w.writerow(["GNAZ_residue","GNAZ_region","RGSZ1_residue_human","RGSZ1_residue_file",
                "RGSZ1_element_Tesmer","contact_type","consistent_type",
                "mindist_rank1","mindist_rank2","mindist_rank3","mean_dist"])
    for r in filt:
        w.writerow([f"{r['arn']}{r['ar']}", gnaz_region(r['ar']),
                    f"{r['brn']}{r['br']+HUMAN}", f"{r['brn']}{r['br']}", rgsz1_element(r['br']),
                    r["consensus"], r["consistent"],
                    f"{r['dists'][0]:.2f}",f"{r['dists'][1]:.2f}",f"{r['dists'][2]:.2f}",f"{r['mean']:.2f}"])
print("wrote",OUT,"|",len(filt),"contacts")
print("%-9s %-10s %-13s %-14s %-13s %s"%("GNAZ","region","RGSZ1(human)","file","Tesmer","type"))
for r in filt:
    print("%-9s %-10s %-13s %-14s %-13s %s"%(f"{r['arn']}{r['ar']}",gnaz_region(r['ar']),
        f"{r['brn']}{r['br']+HUMAN}",f"{r['brn']}{r['br']}",rgsz1_element(r['br']),r['consensus']))
