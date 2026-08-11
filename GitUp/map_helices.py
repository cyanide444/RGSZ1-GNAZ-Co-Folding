#!/usr/bin/env python3
"""Map the 1AGR RGS4 RGS-box secondary structure (deposited HELIX records = Tesmer et al.
1997 assignment) onto RGSZ1's RGS domain via the RGS4->RGSZ1 alignment."""
import numpy as np, re
# 9 helices from 1agr_mod.pdb HELIX records for chain E (RGS4 numbering), N->C = a1..a9
HELICES=[("a1",53,61),("a2",63,68),("a3",70,82),("a4",86,100),("a5",106,118),
         ("a6",131,142),("a7",150,162),("a8",164,170),("a9",172,175)]
# Tesmer 1997 subdomains
SUBDOM={"a1":"terminal","a2":"terminal","a3":"terminal","a8":"terminal","a9":"terminal",
        "a4":"bundle","a5":"bundle","a6":"bundle","a7":"bundle"}
HRGS4="MCKGLAGLPASCLRSAKDMKHRLGFLLQKSDSCEHNSSHNKKDKVVICQRVSQEEVKKWAESLENLISHECGLAAFKAFLKSEYSEENIDFWISCEEYKKIKSPSKLSPKAKKIYNEFISVQATKEVNLDSCTREETSRNMLEPTITCFDEAQKKIFNLMEKDSYRRFLKSRFYLDLVNPSSCGAEKQKGAKSSADCASLVPQCA"
RGSZ1="AQSFDKLMVTPAGRNAFREFLRTEFSEENMLFWMACEELKKEANKNIIEEKARIIYEDYISILSPKEVSLDSRVREVINRNMVEPSQHIFDDAQLQIYTLMHRDSYPRFMNSAVYKDLLQS"
AA1={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE',
 'L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
# ---- BLOSUM62 NW (same as align_1agr.py) ----
_aas="ARNDCQEGHILKMFPSTWYV"
_bl="""4-1-2-2 0-1-1 0-2-1-1-1-1-2-1 1 0-3-2 0
-1 5 0-2-3 1 0-2 0-3-2 2-1-3-2-1-1-3-2-3
-2 0 6 1-3 0 0 0 1-3-3 0-2-3-2 1 0-4-2-3
-2-2 1 6-3 0 2-1-1-3-4-1-3-3-1 0-1-4-3-3
0-3-3-3 9-3-4-3-3-1-1-3-1-2-3-1-1-2-2-1
-1 1 0 0-3 5 2-2 0-3-2 1 0-3-1 0-1-2-1-2
-1 0 0 2-4 2 5-2 0-3-3 1-2-3-1 0-1-3-2-2
0-2 0-1-3-2-2 6-2-4-4-2-3-3-2 0-2-2-3-3
-2 0 1-1-3 0 0-2 8-3-3-1-2-1-2-1-2-2 2-3
-1-3-3-3-1-3-3-4-3 4 2-3 1 0-3-2-1-3-1 3
-1-2-3-4-1-2-3-4-3 2 4-2 2 0-3-2-1-2-1 1
-1 2 0-1-3 1 1-2-1-3-2 5-1-3-1 0-1-3-2-2
-1-1-2-3-1 0-2-3-2 1 2-1 5 0-2-1-1-1-1 1
-2-3-3-3-2-3-3-3-1 0 0-3 0 6-4-2-2 1 3-1
-1-2-2-1-3-1-1-2-2-3-3-1-2-4 7-1-1-4-3-2
1-1 1 0-1 0 0 0-1-2-2 0-1-2-1 4 1-3-2-2
0-1 0-1-1-1-1-2-2-1-1-1-1-2-1 1 5-2-2 0
-3-3-4-4-2-2-3-2-2-3-2-3-1 1-4-3-2 11 2-3
-2-2-2-3-2-1-2-3 2-1-1-2-1 3-3-2-2 2 7-1
0-3-3-3-1-2-2-3-3 3 1-2 1-1-2-2 0-3-1 4"""
rows=[re.findall(r'-?\d+',r) for r in _bl.strip().split("\n")]
B={}
for i,a in enumerate(_aas):
    for j,b in enumerate(_aas): B[(a,b)]=int(rows[i][j])
def sc(a,b): return B.get((a,b),B.get((b,a),-2))
def nw(a,b,gap=-11):
    n,m=len(a),len(b);F=np.zeros((n+1,m+1))
    for i in range(1,n+1):F[i][0]=i*gap
    for j in range(1,m+1):F[0][j]=j*gap
    for i in range(1,n+1):
        for j in range(1,m+1):
            F[i][j]=max(F[i-1][j-1]+sc(a[i-1],b[j-1]),F[i-1][j]+gap,F[i][j-1]+gap)
    i,j=n,m;A2='';B2=''
    while i>0 and j>0:
        if F[i][j]==F[i-1][j-1]+sc(a[i-1],b[j-1]):A2=a[i-1]+A2;B2=b[j-1]+B2;i-=1;j-=1
        elif F[i][j]==F[i-1][j]+gap:A2=a[i-1]+A2;B2='-'+B2;i-=1
        else:A2='-'+A2;B2=b[j-1]+B2;j-=1
    while i>0:A2=a[i-1]+A2;B2='-'+B2;i-=1
    while j>0:A2='-'+A2;B2=b[j-1]+B2;j-=1
    return A2,B2
# align RGS-box core (RGS4 from residue 51) to RGSZ1
A2,B2=nw(HRGS4[50:],RGSZ1)
r2z={}; ra=50; rb=0
for ca,cb in zip(A2,B2):
    if ca!='-':ra+=1
    if cb!='-':rb+=1
    if ca!='-' and cb!='-': r2z[ra]=rb   # RGS4 resnum -> RGSZ1 file resnum
def mapres(r4):
    if r4 in r2z: return r2z[r4], False
    for off in range(1,6):
        if r4-off in r2z: return r2z[r4-off]+off, True
        if r4+off in r2z: return r2z[r4+off]-off, True
    return None, True
def zres(rz): return AA1.get(RGSZ1[rz-1],'?') if rz and 0<rz<=len(RGSZ1) else '---'
def r4res(r4): return AA1.get(HRGS4[r4-1],'?') if 0<r4<=len(HRGS4) else '---'

print("RGS-box secondary structure: 1AGR RGS4  ->  RGSZ1 RGS domain")
print("(RGS4/human-RGS4 numbering; RGSZ1 file numbering 1-121, human = file+261)")
print("%-4s %-9s %-18s %-22s %-9s"%("elem","subdomain","RGS4 range","RGSZ1 file (human)","note"))
print("-"*74)
prev_end=None; prev_name=None
def emit_loop(n1,e4,n2,s4):
    # loop between helix ending e4 (RGS4) and helix starting s4
    l4=(e4+1,s4-1)
    if l4[0]>l4[1]:
        print("%-4s %-9s %-18s %-22s"%(f"{n1}-{n2}","loop","(none)","(direct)")); return
    z1,f1=mapres(l4[0]); z2,f2=mapres(l4[1])
    zr=f"{z1}-{z2} ({z1+261}-{z2+261})" if z1 and z2 else "unmapped"
    flag="~indel" if (f1 or f2) else ""
    print("%-4s %-9s %-18s %-22s %-9s"%(f"{n1}-{n2}","loop",f"{l4[0]}-{l4[1]}",zr,flag))
for name,s4,e4 in HELICES:
    if prev_end is not None: emit_loop(prev_name,prev_end,name,s4)
    zs,fs=mapres(s4); ze,fe=mapres(e4)
    zr=f"{zs}-{ze} ({zs+261}-{ze+261})" if zs and ze else "unmapped"
    flag="~indel" if (fs or fe) else ""
    print("%-4s %-9s %-18s %-22s %-9s"%(name,SUBDOM[name],
          f"{r4res(s4)}{s4}-{r4res(e4)}{e4}", zr, flag))
    prev_end=e4; prev_name=name
print("\nGa-contact loops (Tesmer 1997): a3-a4, a5-a6, a7-a8")
