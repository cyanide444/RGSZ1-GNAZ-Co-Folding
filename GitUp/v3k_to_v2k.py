#!/usr/bin/env python3
"""Faithful V3000 -> V2000 molfile/SDF conversion. Preserves atom order/numbering,
coordinates, bonds, bond wedges (CFG->stereo) and formal charges (via M CHG) exactly."""
import re
SRC=r"C:\Users\richi\Documents\ionizedGTP.sdf"
OUT=r"C:\Users\richi\Documents\ionizedGTP_2.sdf"
lines=open(SRC).read().splitlines()

# header (first 3 lines): title, program, comment
title,prog,comment=lines[0],lines[1],lines[2]

atoms=[]   # (idx, elem, x, y, z, charge)
bonds=[]   # (idx, order, a1, a2, cfg)
mode=None
for l in lines:
    s=l.strip()
    if s.startswith("M  V30 COUNTS"):
        p=s.split(); na=int(p[3]); nb=int(p[4]); chiral=int(p[7]) if len(p)>7 else 0
    elif s=="M  V30 BEGIN ATOM": mode="atom"
    elif s=="M  V30 END ATOM": mode=None
    elif s=="M  V30 BEGIN BOND": mode="bond"
    elif s=="M  V30 END BOND": mode=None
    elif s.startswith("M  V30") and mode=="atom":
        t=s[6:].split()          # idx elem x y z aamap [CHG=..] ...
        idx=int(t[0]); elem=t[1]; x=float(t[2]); y=float(t[3]); z=float(t[4])
        chg=0
        for tok in t[5:]:
            if tok.startswith("CHG="): chg=int(tok.split("=")[1])
        atoms.append((idx,elem,x,y,z,chg))
    elif s.startswith("M  V30") and mode=="bond":
        t=s[6:].split()          # idx order a1 a2 [CFG=..]
        idx=int(t[0]); order=int(t[1]); a1=int(t[2]); a2=int(t[3])
        cfg=0
        for tok in t[4:]:
            if tok.startswith("CFG="): cfg=int(tok.split("=")[1])
        bonds.append((idx,order,a1,a2,cfg))

atoms.sort(key=lambda a:a[0])    # ensure output order == atom index (numbering preserved)
assert [a[0] for a in atoms]==list(range(1,len(atoms)+1)), "atom indices not 1..N sequential!"

CFG2V2K={0:0,1:1,2:4,3:6}        # V3000 CFG -> V2000 bond stereo (up/either/down)
out=[]
out.append(title)
out.append(prog)
out.append(comment)
out.append("%3d%3d  0  0%3d  0  0  0  0  0999 V2000"%(na,nb,chiral))
for idx,elem,x,y,z,chg in atoms:
    out.append("%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0"%(x,y,z,elem))
for idx,order,a1,a2,cfg in bonds:
    out.append("%3d%3d%3d%3d"%(a1,a2,order,CFG2V2K[cfg]))
charged=[(idx,chg) for idx,elem,x,y,z,chg in atoms if chg!=0]
if charged:
    out.append("M  CHG%3d"%len(charged)+"".join("%4d%4d"%(i,c) for i,c in charged))
out.append("M  END")
out.append("$$$$")
open(OUT,"w",newline="\n").write("\n".join(out)+"\n")
print("wrote",OUT)
print("atoms=%d bonds=%d chiral=%d charged atoms=%s"%(na,nb,chiral,[i for i,c in charged]))
