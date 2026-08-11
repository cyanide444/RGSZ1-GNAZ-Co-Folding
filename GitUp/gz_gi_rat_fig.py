#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, re, html, numpy as np
GA=os.path.dirname(os.path.abspath(__file__))
def load(p): return ''.join(l.strip() for l in open(p) if not l.startswith('>'))
gnaz=load(os.path.join(GA,'P19086.fasta')); hgi=load(os.path.join(GA,'P63096.fasta')); rgi=load(os.path.join(GA,'P10824.fasta'))
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
def nw(a,b,gap=-10):
    n,m=len(a),len(b);F=np.zeros((n+1,m+1));P=np.zeros((n+1,m+1),int)
    for i in range(1,n+1):F[i][0]=i*gap;P[i][0]=1
    for j in range(1,m+1):F[0][j]=j*gap;P[0][j]=2
    for i in range(1,n+1):
        for j in range(1,m+1):
            c=[F[i-1][j-1]+sc(a[i-1],b[j-1]),F[i-1][j]+gap,F[i][j-1]+gap]
            k=int(np.argmax(c));F[i][j]=c[k];P[i][j]=k
    i,j=n,m;A='';C=''
    while i>0 or j>0:
        k=P[i][j]
        if i>0 and j>0 and k==0:A=a[i-1]+A;C=b[j-1]+C;i-=1;j-=1
        elif i>0 and k==1:A=a[i-1]+A;C='-'+C;i-=1
        else:A='-'+A;C=b[j-1]+C;j-=1
    return A,C
# two pairwise alignments sharing human GNAI1 (hgi), then merge on the common sequence
Az,Ai=nw(gnaz,hgi)      # GNAZ vs hGNAI1
Bi,Br=nw(hgi,rgi)       # hGNAI1 vs ratGNAI1
def pid2(X,Y):
    idn=al=0
    for x,y in zip(X,Y):
        if x=='-' or y=='-': continue
        al+=1; idn+=(x==y)
    return 100.0*idn/al, idn, al
PID_gz,idn_gz,al_gz=pid2(Az,Ai)
PID_rat,idn_rat,al_rat=pid2(Bi,Br)
# merge: walk A (GNAZ/hgi) and B (hgi/rat), synchronising on hgi residues
p=q=0; mgz=[]; mhi=[]; mrat=[]
while p<len(Ai) or q<len(Bi):
    ap_gap = p<len(Ai) and Ai[p]=='-'
    bq_gap = q<len(Bi) and Bi[q]=='-'
    if ap_gap:                        # GNAZ insertion vs hgi
        mgz.append(Az[p]); mhi.append('-'); mrat.append('-'); p+=1
    elif bq_gap:                      # rat insertion vs hgi
        mgz.append('-'); mhi.append('-'); mrat.append(Br[q]); q+=1
    else:                             # both consume an hgi residue
        mgz.append(Az[p]); mhi.append(Bi[q]); mrat.append(Br[q]); p+=1; q+=1
# per-column numbering
gz_n=hi_n=rt_n=0; col=[]
for x,y,z in zip(mgz,mhi,mrat):
    gz=hi=rt=None
    if x!='-': gz_n+=1; gz=gz_n
    if y!='-': hi_n+=1; hi=hi_n
    if z!='-': rt_n+=1; rt=rt_n
    col.append((x,y,z,gz))
SW={'Switch I':(175,188,'#f6c96b'),'Switch II':(200,220,'#8fca8f'),'Switch III':(228,240,'#9dc0ee')}
def swof(gz):
    if gz is None: return None
    for nm,(lo,hi,cc) in SW.items():
        if lo<=gz<=hi: return cc
    return None
SEQIDX={'GNAZ':0,'hGI':1,'rGI':2}
per=60; nblk=(len(col)+per-1)//per
def cells(which,c0,c1):
    out=[]
    for k in range(c0,c1):
        ch=col[k][SEQIDX[which]]; cc=swof(col[k][3])
        st='background:%s;'%cc if cc else ''
        if ch=='-': st+='color:#bbb;'
        out.append('<span style="%s">%s</span>'%(st,html.escape(ch)))
    return ''.join(out)
def consline(c0,c1):
    out=[]
    for k in range(c0,c1):
        a,b,c,_=col[k]
        if '-' in (a,b,c): ch='&nbsp;'
        elif a==b==c: ch='*'
        elif sc(a,b)>0 and sc(a,c)>0 and sc(b,c)>0: ch=':'
        else: ch='&nbsp;'
        out.append(ch)
    return ''.join(out)
def num(which,c0,c1,last=False):
    idx=SEQIDX[which]; vals=[]
    for k in range(c0,c1):
        x,y,z,_=col[k]; v=(x,y,z)[idx]
        # recover running number
    # simpler: track numbers directly
    return ''
# track running numbers per sequence for block endpoints
runs={w:[] for w in SEQIDX}
cn={'GNAZ':0,'hGI':0,'rGI':0}
for k in range(len(col)):
    a,b,c,_=col[k]
    for w,ch in (('GNAZ',a),('hGI',b),('rGI',c)):
        if ch!='-': cn[w]+=1
        runs[w].append(cn[w] if ch!='-' else None)
def first(w,c0,c1):
    for k in range(c0,c1):
        if runs[w][k] is not None: return runs[w][k]
    return ''
def lastn(w,c0,c1):
    v=''
    for k in range(c0,c1):
        if runs[w][k] is not None: v=runs[w][k]
    return v
LAB={'GNAZ':'GNAZ (human)','hGI':'GNAI1 (human)','rGI':'GNAI1 (rat)'}
blocks=[]
for bk in range(nblk):
    c0=bk*per; c1=min(c0+per,len(col)); r=[]
    for w in ('GNAZ','hGI','rGI'):
        r.append('<div class="row"><span class="lab">%s</span><span class="num">%4s</span><span class="sq">%s</span><span class="num">&nbsp;%-4s</span></div>'
                 %(LAB[w],first(w,c0,c1),cells(w,c0,c1),lastn(w,c0,c1)))
    r.append('<div class="row"><span class="lab"></span><span class="num"></span><span class="mid">%s</span></div>'%consline(c0,c1))
    blocks.append('<div class="blk">'+''.join(r)+'</div>')
leg=' '.join('<span class="lg"><span class="sw" style="background:%s"></span>%s (GNAZ %d&ndash;%d)</span>'%(cc,nm,lo,hi) for nm,(lo,hi,cc) in SW.items())
BODY='''<div class="wrap">
<h1>Human G&alpha;z (GNAZ) aligned with human and rat G&alpha;i1 (GNAI1)</h1>
<p class="sub">Multiple alignment of full-length UniProt sequences: human GNAZ (P19086), human GNAI1 (P63096) and
rat GNAI1 (P10824). <b>*</b> = identical in all three sequences, <b>:</b> = conservative in all three
(all pairwise BLOSUM62 scores positive). The three switch regions (GNAZ numbering) are shaded.
GNAZ is %.1f%% identical to human GNAI1; human and rat GNAI1 are %.1f%% identical.</p>
<div class="legend">%s</div>
<div class="mono">%s</div>
<div class="cap">Alignment built from two pairwise global alignments (Needleman&ndash;Wunsch, BLOSUM62), each sharing
the human GNAI1 sequence, then merged on the common sequence. Switch&nbsp;I 175&ndash;188 (&gamma;-phosphate/Mg<sup>2+</sup>
sensor Thr182); Switch&nbsp;II 200&ndash;220 (catalytic Gln205, Arg-finger region); Switch&nbsp;III 228&ndash;240.
Rat GNAI1 (P10824) is the sequence used in the RGS4&middot;Gi1 crystal structure 1AGR.</div>
</div>'''%(PID_gz,PID_rat,leg,''.join(blocks))
STYLE='''*{box-sizing:border-box}body{margin:0;background:#fff}
.wrap{padding:24px 30px;font-family:'Segoe UI',system-ui,sans-serif;color:#111;display:inline-block}
h1{font-size:19px;margin:0 0 6px}.sub{font-size:12.5px;color:#444;margin:0 0 10px;max-width:1000px;line-height:1.45}
.legend{margin:4px 0 16px;font-size:12px;color:#333}
.lg{margin-right:22px;white-space:nowrap}.sw{display:inline-block;width:13px;height:13px;border:1px solid #999;vertical-align:-2px;margin-right:5px}
.mono{font-family:'Consolas','Courier New',monospace;font-size:15px;line-height:1.15;white-space:pre}
.blk{margin-bottom:16px}.row{white-space:pre}
.lab{color:#333;font-weight:700;display:inline-block;width:120px}
.num{color:#999;display:inline-block;width:40px;text-align:right}
.sq span{padding:1px 0;letter-spacing:2px}.mid{color:#888;letter-spacing:2px}
.cap{margin-top:8px;font-size:11px;color:#666;max-width:1000px;line-height:1.5;border-top:1px solid #e3e3dd;padding-top:8px}'''
open(os.path.join(GA,'_gzgirat.html'),'w',encoding='utf-8').write('<!doctype html><html><head><meta charset="utf-8"><style>%s</style></head><body>%s</body></html>'%(STYLE,BODY))
print('GNAZ vs hGNAI1 = %.1f%% (%d/%d); hGNAI1 vs rGNAI1 = %.1f%% (%d/%d); cols=%d'%(PID_gz,idn_gz,al_gz,PID_rat,idn_rat,al_rat,len(col)))
