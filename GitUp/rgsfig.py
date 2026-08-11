#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re, os, html, numpy as np
D=os.path.dirname(os.path.abspath(__file__))
DOM={'Q08116':(85,200),'P41220':(83,199),'P49796':(1073,1198),'P49798':(62,178),'O15539':(64,180),
 'P57771':(56,171),'O14921':(34,150),'O15492':(65,181),'Q9NS28':(86,202),'Q2M5E4':(21,137),
 'P49758':(336,441),'P49802':(333,448),'O75916':(298,413),'O94810':(299,414),
 'Q9UGC6':(84,200),'P49795':(90,206),'O76081':(262,378)}
INFO=[('RGS1','Q08116','R4'),('RGS2','P41220','R4'),('RGS3','P49796','R4'),('RGS4','P49798','R4'),
 ('RGS5','O15539','R4'),('RGS8','P57771','R4'),('RGS13','O14921','R4'),('RGS16','O15492','R4'),
 ('RGS18','Q9NS28','R4'),('RGS21','Q2M5E4','R4'),('RGS6','P49758','R7'),('RGS7','P49802','R7'),
 ('RGS9','O75916','R7'),('RGS11','O94810','R7'),('RGS17','Q9UGC6','RZ'),('RGS19','P49795','RZ'),('RGS20','O76081','RZ')]
seq={a:''.join(l.strip() for l in open(os.path.join(D,a+'.fasta')) if not l.startswith('>')) for _,a,_ in INFO}
boxs={a:seq[a][DOM[a][0]-1:DOM[a][1]] for _,a,_ in INFO}
# BLOSUM62 NW
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
def nw(a,b,gap=-9):
    n,m=len(a),len(b);F=np.zeros((n+1,m+1))
    for i in range(1,n+1):F[i][0]=i*gap
    for j in range(1,m+1):F[0][j]=j*gap
    for i in range(1,n+1):
        for j in range(1,m+1):
            F[i][j]=max(F[i-1][j-1]+sc(a[i-1],b[j-1]),F[i-1][j]+gap,F[i][j-1]+gap)
    i,j=n,m;A='';C=''
    while i>0 and j>0:
        if F[i][j]==F[i-1][j-1]+sc(a[i-1],b[j-1]):A=a[i-1]+A;C=b[j-1]+C;i-=1;j-=1
        elif F[i][j]==F[i-1][j]+gap:A=a[i-1]+A;C='-'+C;i-=1
        else:A='-'+A;C=b[j-1]+C;j-=1
    while i>0:A=a[i-1]+A;C='-'+C;i-=1
    while j>0:A='-'+A;C=b[j-1]+C;j-=1
    return A,C
ref=boxs['P49798']; L=len(ref); CAT=128-62      # catalytic column index (RGS4 N128)
def project(sb):
    a,b=nw(ref,sb); out=[]
    for x,y in zip(a,b):
        if x!='-': out.append(y)
    return out
msa={a:project(boxs[a]) for _,a,_ in INFO}
# RGS4 helices -> box-index ranges (color the different helices)
HEL=[(1,6),(8,20),(24,38),(44,56),(69,80),(88,100),(102,108),(110,113)]
HCOL=['#f4a6b8','#f9c48a','#f4e58c','#bfe3a0','#99e0c8','#9fd4ea','#aeb8ee','#cdb4e8']
def hcol(idx):
    for (a,b),c in zip(HEL,HCOL):
        if a<=idx<=b: return c
    return None
FAMBG={'R4':'#f7f7f4','R7':'#f2f6fa','RZ':'#fbf1f4'}
FAMLBL={'R4':'R4 subfamily','R7':'R7 subfamily','RZ':'RZ subfamily (Gz / Gi-selective)'}

per=60; blocks=[(s,min(s+per,L)) for s in range(0,L,per)]
def render():
    out=[]
    for bi,(c0,c1) in enumerate(blocks):
        # marker row: N128 label above catalytic column if in this block
        if c0<=CAT<c1:
            pos=CAT-c0; mk=['&nbsp;']*(c1-c0)
            lab='N128'; start=min(pos,(c1-c0)-len(lab))
            for k,ch in enumerate(lab): mk[start+k]=ch
            out.append('<div class="mk">%s%s</div>'%('&nbsp;'*10,''.join(mk)))
        prevfam=None
        for name,acc,fam in INFO:
            if fam!=prevfam:
                out.append('<div class="fam">%s</div>'%FAMLBL[fam]); prevfam=fam
            cells=['<span class="lab" style="background:%s">%-3s %-6s</span>'%(FAMBG[fam],fam,name)]
            for idx in range(c0,c1):
                ch=msa[acc][idx]; st=''
                c=hcol(idx)
                if c: st+='background:%s;'%c
                cat=(idx==CAT)
                if cat:
                    st+='outline:2.5px solid #d00;outline-offset:-1px;font-weight:700;'
                    if ch=='S': st+='color:#d00;'
                    elif ch=='N': st+='color:#111;'
                cells.append('<span style="%s">%s</span>'%(st, html.escape(ch)))
            out.append('<div class="seq" style="background:%s">%s</div>'%(FAMBG[fam],''.join(cells)))
        out.append('<div class="sp"></div>')
    return ''.join(out)

BODY=f'''<div class="wrap">
<h1>The invariant RGS catalytic asparagine (RGS4 Asn128) is a serine in the RZ subfamily</h1>
<p class="sub">RGS domains of all human R4, R7 and RZ subfamily members (UniProt), aligned to RGS4. The eight
RGS-box helices are highlighted in color; the catalytic column (RGS4 Asn128) is boxed in red.
<b>Asn is invariant in R4 (10/10) and R7 (4/4); all three RZ members (RGS17, RGS19, RGS20) carry a serine.</b></p>
<div class="mono">{render()}</div>
<div class="cap">Aligned by reference projection onto the RGS4 RGS box (BLOSUM62 Needleman&ndash;Wunsch); columns
correspond to RGS4 residues 62&ndash;178. Colored blocks mark the eight &alpha;-helices of the RGS fold.
Red&nbsp;box = the catalytic position; the RZ serines are shown in red.</div>
</div>'''
STYLE='''*{box-sizing:border-box}body{margin:0;background:#fff}
.wrap{padding:24px 28px;font-family:'Segoe UI',system-ui,sans-serif;color:#111;display:inline-block}
h1{font-size:18px;margin:0 0 5px}.sub{font-size:12.5px;color:#444;margin:0 0 14px;max-width:1050px;line-height:1.45}
.mono{font-family:'Consolas','Courier New',monospace;font-size:14px;line-height:1.28;white-space:pre}
.seq{letter-spacing:1px}.mono span{padding:1px 0}
.lab{color:#222;font-weight:600;letter-spacing:0}
.fam{font-family:'Segoe UI',sans-serif;font-weight:700;font-size:11.5px;color:#555;margin:7px 0 1px}
.mk{color:#d00;font-weight:700;letter-spacing:1px;height:16px}
.sp{height:12px}
.cap{margin-top:12px;font-size:11px;color:#666;max-width:1050px;line-height:1.5;border-top:1px solid #e3e3dd;padding-top:8px}'''
open(os.path.join(D,'_rgsfig.html'),'w',encoding='utf-8').write('<!doctype html><html><head><meta charset="utf-8"><style>%s</style></head><body>%s</body></html>'%(STYLE,BODY))
print('catalytic column residues:', {n:msa[a][CAT] for n,a,f in INFO})
