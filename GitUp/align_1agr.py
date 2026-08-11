#!/usr/bin/env python3
import json, os, numpy as np
D=r"C:\Users\richi\Documents\1AGR"
J=json.load(open(os.path.join(D,"_iface.json")))
AA1={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
 'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
# reference sequences
HGI ="MGCTLSAEDKAAVERSKMIDRNLREDGEKAAREVKLLLLGAGESGKSTIVKQMKIIHEAGYSEEECKQYKAVVYSNTIQSIIAIIRAMGRLKIDFGDSARADDARQLFVLAGAAEEGFMTAELAGVIKRLWKDSGVQACFNRSREYQLNDSAAYYLNDLDRIAQPNYIPTQQDVLRTRVKTTGIVETHFTFKDLHFKMFDVGGQRSERKKWIHCFEGVTAIIFCVALSDYDLVLAEDEEEMNRMHESMKLFDSICNNKWFTDTSIILFLNKKDLFEEKIKKSPLTICYPEYAGSNTYEEAAAYIQCQFEDLNKRKDTKEIYTHFTCATDTKNVQFVFDAVTDVIIKNNLKDCGLF"
HRGS4="MCKGLAGLPASCLRSAKDMKHRLGFLLQKSDSCEHNSSHNKKDKVVICQRVSQEEVKKWAESLENLISHECGLAAFKAFLKSEYSEENIDFWISCEEYKKIKSPSKLSPKAKKIYNEFISVQATKEVNLDSCTREETSRNMLEPTITCFDEAQKKIFNLMEKDSYRRFLKSRFYLDLVNPSSCGAEKQKGAKSSADCASLVPQCA"
GNAZ ="MGCRQSSEEKEAARRSRRIDRHLRSESQRQRREIKLLLLGTSNSGKSTIVKQMKIIHSGGFNLEACKEYKPLIIYNAIDSLTRIIRALAALRIDFHNPDRAYDAVQLFALTGPAESKGEITPELLGVMRRLWADPGAQACFSRSSEYHLEDNAAYYLNDLERIAAADYIPTVEDILRSRDMTTGIVENKFTFKELTFKMVDVGGQRSERKKWIHCFEGVTAIIFCVELSGYDLKLYEDNQTSRMAESLRLFDSICNNNWFINTSLILFLNKKDLLAEKIRRIPLTICFPEYKGQNTYEEAAVYIQRQFEDLNRNKETKEIYSHFTCATDTSNIQFVFDAVTDVIIQNNLKYIGLC"
RGSZ1="AQSFDKLMVTPAGRNAFREFLRTEFSEENMLFWMACEELKKEANKNIIEEKARIIYEDYISILSPKEVSLDSRVREVINRNMVEPSQHIFDDAQLQIYTLMHRDSYPRFMNSAVYKDLLQS"
# 1AGR observed sequences (from resolved residues)
def obs(seqd):
    s="".join(AA1.get(seqd[k],'X') for k in sorted(seqd,key=int))
    nums=sorted(int(k) for k in seqd); return s,nums
A_GI,_=obs({int(k):v for k,v in J['gi_seq'].items()})
A_RGS,_=obs({int(k):v for k,v in J['rgs_seq'].items()})
gi_seq={int(k):v for k,v in J['gi_seq'].items()}; rgs_seq={int(k):v for k,v in J['rgs_seq'].items()}
# ---- BLOSUM62 NW ----
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
import re
rows=[re.findall(r'-?\d+',r) for r in _bl.strip().split("\n")]
B={}
for i,a in enumerate(_aas):
    for j,b in enumerate(_aas): B[(a,b)]=int(rows[i][j])
def sc(a,b): return B.get((a,b),B.get((b,a),-2))
def nw(a,b,gap=-8):
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
def mapnums(a,b,a0=1,b0=1,gap=-8):
    # returns dict a_resnum -> (b_char, b_resnum)
    A2,B2=nw(a,b,gap); mp={};ra=a0-1;rb=b0-1
    for ca,cb in zip(A2,B2):
        if ca!='-':ra+=1
        if cb!='-':rb+=1
        if ca!='-' and cb!='-': mp[ra]=(cb,rb)
    return mp,A2,B2
# human Gai1 -> GNAZ ; human RGS4 -> RGSZ1
gi2gnaz,alnGIa,alnGIb=mapnums(HGI,GNAZ)
# RGS4 has a ~50-res N-terminal extension absent from the RGSZ1 1-121 fragment; align the
# RGS-box core only (from residue 51) with a stiffer gap penalty to avoid spurious terminal indels.
rgs2rz,alnRGa,alnRGb=mapnums(HRGS4[50:],RGSZ1,a0=51,gap=-11)
def blocks(a,b,w=70):
    out=[];ia=0
    for k in range(0,len(a),w):
        sa=a[k:k+w];sb=b[k:k+w]
        mk="".join('|' if x==y and x!='-' else (':' if x!='-' and y!='-' and sc(x,y)>0 else ' ') for x,y in zip(sa,sb))
        out+=[sa,mk,sb,""]
    return "\n".join(out)

# ---- restraint sets (refined) ----
gi_active=[116,180,182,184,207,210,235,236]; gi_passive=[185,209,213,237]
rgs_active=[83,84,87,88,128,131,134,161,163,166]; rgs_passive=[124,159,167]

def cons_row(N,seqd,human,mp,partner):
    a1=AA1.get(seqd.get(N,'---'),'-')          # 1AGR AA
    hu=human[N-1] if 0<N<=len(human) else '-'   # human AA (same numbering)
    c1=(a1==hu)
    if N in mp: pc,pn=mp[N]; cg=(pc==hu)
    else: pc,pn,cg='-','-',False
    return a1,hu,c1,pc,pn,cg

# GNAZ table
print("="*78)
print("GNAZ active/passive restraints (mapped from Gai1 via 1AGR)")
print("%-6s %-5s | %-8s | %-14s | %-10s"%("class","Gai1","1AGR/hum","-> GNAZ res","conserved?"))
gnaz_out=[]
for cls,lst in [("ACTIVE",gi_active),("PASSIVE",gi_passive)]:
    for N in lst:
        a1,hu,c1,pc,pn,cg=cons_row(N,gi_seq,HGI,gi2gnaz,"GNAZ")
        gz="%s%s"%(pc,pn) if pn!='-' else "(gap)"
        print("%-6s %-5d | 1AGR=%s hum=%s | GNAZ %-8s | %s"%(cls,N,a1,hu,gz,"conserved" if cg else ("similar" if pn!='-' else "no-map")))
        gnaz_out.append((cls,N,a1,hu,pc,pn,cg))
print("\n"+"="*78)
print("RGSZ1 active/passive restraints (mapped from RGS4 via 1AGR)")
for cls,lst in [("ACTIVE",rgs_active),("PASSIVE",rgs_passive)]:
    for N in lst:
        a1,hu,c1,pc,pn,cg=cons_row(N,rgs_seq,HRGS4,rgs2rz,"RGSZ1")
        rz="%s%s"%(pc,pn) if pn!='-' else "(gap)"
        print("%-6s RGS4 %-4d | 1AGR=%s hum=%s | RGSZ1 %-8s | %s"%(cls,N,a1,hu,rz,"conserved" if cg else ("similar" if pn!='-' else "no-map")))

# ---- write alignment file ----
with open(os.path.join(D,"1AGR_alignments.txt"),"w",encoding="utf-8") as f:
    f.write("SEQUENCE ALIGNMENTS FOR 1AGR-BASED RESTRAINT MAPPING\n"+"="*70+"\n")
    f.write("(| identical, : similar/positive BLOSUM62)\n\n")
    f.write("### 1. 1AGR Gai1 (chain A, resolved) vs human Gai1 (P63096)\n")
    a2,b2=nw(A_GI,HGI); f.write(blocks(a2,b2)+"\n")
    f.write("### 2. human Gai1 (P63096) vs GNAZ (P19086)\n"+blocks(alnGIa,alnGIb)+"\n")
    f.write("### 3. 1AGR RGS4 (chain E, resolved) vs human RGS4 (P49798)\n")
    a2,b2=nw(A_RGS,HRGS4); f.write(blocks(a2,b2)+"\n")
    f.write("### 4. human RGS4 (P49798) vs RGSZ1 RGS domain (docked 1-121)\n"+blocks(alnRGa,alnRGb)+"\n")
print("\nWrote 1AGR_alignments.txt")

# ---- write conservation file ----
with open(os.path.join(D,"1AGR_conservation_GNAZ_RGSZ1.txt"),"w",encoding="utf-8") as f:
    f.write("CONSERVATION OF PROPOSED 1AGR RESTRAINT RESIDUES IN GNAZ / RGSZ1\n"+"="*70+"\n\n")
    f.write("Gai1 -> GNAZ (Gai1 numbering; 1AGR & human Gai1 share numbering)\n")
    f.write("%-6s %-8s %-8s %-14s %-12s\n"%("class","Gai1","1AGR=hum?","-> GNAZ","conserved in GNAZ?"))
    for cls,lst in [("active",gi_active),("passive",gi_passive)]:
        for N in lst:
            a1,hu,c1,pc,pn,cg=cons_row(N,gi_seq,HGI,gi2gnaz,"GNAZ")
            f.write("%-6s %s%-6d %-8s GNAZ %s%-8s %s\n"%(cls,a1,N,("yes" if c1 else "NO(%s)"%hu),pc,pn,("YES (%s=%s)"%(pc,hu) if cg else ("no: GNAZ=%s vs %s"%(pc,hu) if pn!='-' else "no map"))))
    f.write("\nRGS4 -> RGSZ1 (RGS4 numbering)\n")
    f.write("%-6s %-8s %-8s %-14s %-12s\n"%("class","RGS4","1AGR=hum?","-> RGSZ1","conserved in RGSZ1?"))
    for cls,lst in [("active",rgs_active),("passive",rgs_passive)]:
        for N in lst:
            a1,hu,c1,pc,pn,cg=cons_row(N,rgs_seq,HRGS4,rgs2rz,"RGSZ1")
            f.write("%-6s %s%-6d %-8s RGSZ1 %s%-8s %s\n"%(cls,a1,N,("yes" if c1 else "NO(%s)"%hu),pc,pn,("YES (%s=%s)"%(pc,hu) if cg else ("no: RGSZ1=%s vs %s"%(pc,hu) if pn!='-' else "no map"))))
print("Wrote 1AGR_conservation_GNAZ_RGSZ1.txt")
