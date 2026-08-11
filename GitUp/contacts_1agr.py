#!/usr/bin/env python3
"""Compute the Gai1(chain A)-RGS4(chain E) interface in 1agr_mod.pdb, assess the
user-proposed restraint residues, and assign active/passive."""
import math, os
F=r"C:\Users\richi\Documents\1AGR\1agr_mod.pdb"
AA1={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
 'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
HYD={'ALA','VAL','LEU','ILE','MET','PHE','TRP','PRO','CYS'}
DET=5.0; SB=4.0; HB=3.5
prop_gi=sorted(set([180,182,210,236,237,116,235,213,185,209,236]))
prop_rgs=sorted(set([84,87,128,163,124,159,134,131,166,167,83,179,161]))
def cat(rn,a):return (rn=='LYS' and a=='NZ') or (rn=='ARG' and a in('NE','NH1','NH2'))
def ani(rn,a):return (rn=='ASP' and a in('OD1','OD2')) or (rn=='GLU' and a in('OE1','OE2'))
def parse(ch):
    res={}
    for l in open(F):
        if l[:6]!='ATOM  ' or l[21]!=ch: continue
        rn=l[17:20].strip(); ri=int(l[22:26]); nm=l[12:16].strip(); el=(l[76:78].strip() or nm[0])
        if el=='H': continue
        res.setdefault(ri,[rn,[]])[1].append((nm,el,(float(l[30:38]),float(l[38:46]),float(l[46:54]))))
    return res
A=parse('A')  # Gai1
E=parse('E')  # RGS4
def dist(p,q):return math.dist(p,q)
def classify(rn1,at1,rn2,at2):
    dm=9e9;salt=hb=hyd=False;ca=None
    for n1,e1,x1 in at1:
        for n2,e2,x2 in at2:
            d=dist(x1,x2)
            if d>DET: continue
            if d<dm: dm=d;ca=(n1,n2)
            if d<=SB and ((cat(rn1,n1) and ani(rn2,n2)) or (ani(rn1,n1) and cat(rn2,n2))): salt=True
            if d<=HB and e1 in('N','O') and e2 in('N','O'): hb=True
            if e1=='C' and e2=='C' and rn1 in HYD and rn2 in HYD and d<=4.5: hyd=True
    if dm>DET: return None
    return dm,ca,('Salt bridge' if salt else 'H-bond' if hb else 'Hydrophobic' if hyd else 'van der Waals')
# per Gai1 residue -> best RGS4 contact
gi_hits={}; rgs_hits={}
for gi in A:
    best=None
    for ei in E:
        r=classify(A[gi][0],A[gi][1],E[ei][0],E[ei][1])
        if r and (best is None or r[0]<best[0]): best=(r[0],ei,r[2])
        if r:
            gi_hits.setdefault(gi,[]).append((ei,r[2],r[0]))
            rgs_hits.setdefault(ei,[]).append((gi,r[2],r[0]))
gi_iface=sorted(gi_hits); rgs_iface=sorted(rgs_hits)
def besttype(hits):
    types=[t for _,t,_ in hits]
    return 'Salt bridge' if 'Salt bridge' in types else 'H-bond' if 'H-bond' in types else 'Hydrophobic' if 'Hydrophobic' in types else 'van der Waals'
def mind(hits): return min(d for _,_,d in hits)
print("=== Gai1 (chain A) interface residues (<=%.1f A to RGS4) ==="%DET)
for gi in gi_iface:
    print("  %s%-4d  %-13s min %.2f  <-in proposal:%s"%(AA1.get(A[gi][0],A[gi][0]),gi,besttype(gi_hits[gi]),mind(gi_hits[gi]),gi in prop_gi))
print("\n=== RGS4 (chain E) interface residues ===")
for ei in rgs_iface:
    print("  %s%-4d  %-13s min %.2f  <-in proposal:%s"%(AA1.get(E[ei][0],E[ei][0]),ei,besttype(rgs_hits[ei]),mind(rgs_hits[ei]),ei in prop_rgs))
print("\n=== assessment ===")
print("proposed Gai1 NOT at interface:", [g for g in prop_gi if g not in gi_hits])
print("interface Gai1 NOT in proposal :", [g for g in gi_iface if g not in prop_gi])
print("proposed RGS4 NOT at interface:", [r for r in prop_rgs if r not in rgs_hits])
print("interface RGS4 NOT in proposal :", [r for r in rgs_iface if r not in prop_rgs])
# active = has H-bond or salt bridge; passive = vdW/hydrophobic only (among interface residues)
def split(iface,hits):
    act=[r for r in iface if besttype(hits[r]) in ('Salt bridge','H-bond')]
    pas=[r for r in iface if r not in act]
    return act,pas
ga,gp=split(gi_iface,gi_hits); ra,rp=split(rgs_iface,rgs_hits)
print("\nGai1 active :", ga)
print("Gai1 passive:", gp)
print("RGS4 active :", ra)
print("RGS4 passive:", rp)
# export residue lists for the alignment script
import json
open(r"C:\Users\richi\Documents\1AGR\_iface.json","w").write(json.dumps(dict(
  gi_iface=gi_iface,rgs_iface=rgs_iface,gi_active=ga,gi_passive=gp,rgs_active=ra,rgs_passive=rp,
  prop_gi=prop_gi,prop_rgs=prop_rgs,
  gi_seq={i:A[i][0] for i in A}, rgs_seq={i:E[i][0] for i in E})))
