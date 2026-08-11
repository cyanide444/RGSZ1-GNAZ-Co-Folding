#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, re, sys, numpy as np
GA=os.path.dirname(os.path.abspath(__file__))
RG=os.path.join(os.path.dirname(GA),'rgsfam')
def load(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>'))
# BLOSUM62 + NW (from rgsfig.py)
_aas='ARNDCQEGHILKMFPSTWYV'
_bl='''4-1-2-2 0-1-1 0-2-1-1-1-1-2-1 1 0-3-2 0
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
0-3-3-3-1-2-2-3-3 3 1-2 1-1-2-2 0-3-1 4'''
rows=[re.findall(r'-?\d+',r) for r in _bl.strip().split(chr(10))];B={}
for i,x in enumerate(_aas):
    for j,y in enumerate(_aas): B[(x,y)]=int(rows[i][j])
sc=lambda a,b:B.get((a,b),B.get((b,a),-2))
def nw(a,b,gap=-10,ext=-0.5):
    n,m=len(a),len(b);F=np.zeros((n+1,m+1))
    for i in range(1,n+1):F[i][0]=gap+(i-1)*ext
    for j in range(1,m+1):F[0][j]=gap+(j-1)*ext
    P=np.zeros((n+1,m+1),int)
    for i in range(1,n+1):
        for j in range(1,m+1):
            c=[F[i-1][j-1]+sc(a[i-1],b[j-1]),F[i-1][j]+gap,F[i][j-1]+gap]
            k=int(np.argmax(c));F[i][j]=c[k];P[i][j]=k
    i,j=n,m;A='';C=''
    while i>0 and j>0:
        k=P[i][j]
        if k==0:A=a[i-1]+A;C=b[j-1]+C;i-=1;j-=1
        elif k==1:A=a[i-1]+A;C='-'+C;i-=1
        else:A='-'+A;C=b[j-1]+C;j-=1
    while i>0:A=a[i-1]+A;C='-'+C;i-=1
    while j>0:A='-'+A;C=b[j-1]+C;j-=1
    return A,C
def pid(a,b):
    A,C=nw(a,b); idn=al=0
    for x,y in zip(A,C):
        if x=='-' or y=='-': continue
        al+=1; idn+=(x==y)
    return 100.0*idn/al

# ---------------- Galpha ----------------
GAL=[('GNAZ','P19086'),('GNAI1','P63096'),('GNAI2','P04899'),('GNAI3','P08754'),('GNAO','P09471')]
gs={n:load(os.path.join(GA,a+'.fasta')) for n,a in GAL}
names=[n for n,_ in GAL]
M=np.zeros((5,5))
for i in range(5):
    for j in range(5):
        M[i,j]=100.0 if i==j else pid(gs[names[i]],gs[names[j]])
print('=== Galpha pairwise %identity (global NW, BLOSUM62) ===')
print('%-7s'%'', ' '.join('%7s'%n for n in names))
for i,n in enumerate(names):
    print('%-7s'%n, ' '.join('%7.1f'%M[i,j] for j in range(5)))
print()

# ---------------- RGS boxes ----------------
DOM={'Q08116':(85,200),'P41220':(83,199),'P49796':(1073,1198),'P49798':(62,178),'O15539':(64,180),
 'P57771':(56,171),'O14921':(34,150),'O15492':(65,181),'Q9NS28':(86,202),'Q2M5E4':(21,137),
 'P49758':(336,441),'P49802':(333,448),'O75916':(298,413),'O94810':(299,414),
 'Q9UGC6':(84,200),'P49795':(90,206),'O76081':(262,378)}
INFO=[('RGS1','Q08116','R4'),('RGS2','P41220','R4'),('RGS3','P49796','R4'),('RGS4','P49798','R4'),
 ('RGS5','O15539','R4'),('RGS8','P57771','R4'),('RGS13','O14921','R4'),('RGS16','O15492','R4'),
 ('RGS18','Q9NS28','R4'),('RGS21','Q2M5E4','R4'),('RGS6','P49758','R7'),('RGS7','P49802','R7'),
 ('RGS9','O75916','R7'),('RGS11','O94810','R7'),('RGS17','Q9UGC6','RZ'),('RGS19','P49795','RZ'),('RGS20','O76081','RZ')]
box={}
for nm,a,f in INFO:
    s=load(os.path.join(RG,a+'.fasta')); box[nm]=s[DOM[a][0]-1:DOM[a][1]]
rn=[nm for nm,_,_ in INFO]; fam={nm:f for nm,_,f in INFO}
K=len(rn); R=np.zeros((K,K))
for i in range(K):
    for j in range(K):
        R[i,j]=100.0 if i==j else pid(box[rn[i]],box[rn[j]])
np.save(os.path.join(GA,'rgs_idmat.npy'),R)
np.save(os.path.join(GA,'gal_idmat.npy'),M)
import json; json.dump({'rn':rn,'fam':[fam[x] for x in rn],'gnames':names}, open(os.path.join(GA,'labels.json'),'w'))
# family summaries
def block(fa,fb):
    vals=[R[i,j] for i in range(K) for j in range(K) if i!=j and fam[rn[i]]==fa and fam[rn[j]]==fb and (fa!=fb or i<j)]
    return np.mean(vals),min(vals),max(vals),len(vals)
print('=== RGS-box pairwise %identity summary ===')
print('%-14s %6s %6s %6s %5s'%('comparison','mean','min','max','npair'))
for fa,fb in [('R4','R4'),('R7','R7'),('RZ','RZ'),('R4','R7'),('R4','RZ'),('R7','RZ')]:
    mu,mn,mx,np_=block(fa,fb)
    print('%-14s %6.1f %6.1f %6.1f %5d'%('%s vs %s'%(fa,fb),mu,mn,mx,np_))
allv=[R[i,j] for i in range(K) for j in range(K) if i<j]
print('%-14s %6.1f %6.1f %6.1f %5d'%('ALL pairs',np.mean(allv),min(allv),max(allv),len(allv)))
# RGS20 (RGSZ1) row
z=rn.index('RGS20')
print('\nRGS20 (RGSZ1) vs others:')
order=sorted([x for x in range(K) if x!=z],key=lambda j:-R[z,j])
for j in order: print('   %-6s (%s) %5.1f'%(rn[j],fam[rn[j]],R[z,j]))
