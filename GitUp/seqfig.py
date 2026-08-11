#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, re, html, numpy as np
SCR=os.path.dirname(os.path.abspath(__file__))
def load(acc):
    return ''.join(l.strip() for l in open(os.path.join(SCR,acc+'.fasta')) if not l.startswith('>'))
GNAZ=load('P19086'); RGS=load('O76081'); GNAI1=load('P63096')

# ---- colors (highlighter-style) ----
C_SW1='#c3d0f5'; C_SW2='#ffd9a8'; C_SW3='#e6c2ef'; C_HELIX='#bfe6c8'; C_DOM='#efefe9'; C_HD='#e8e0cc'
INK='#111'; SEC='#555'; NUMC='#999'

# ---- annotations ----
# GNAZ switches (residue -> class)
GNAZ_ANN={}
for a,b,c in [(63,174,'hd'),(175,188,'sw1'),(200,220,'sw2'),(228,240,'sw3')]:
    for i in range(a,b+1): GNAZ_ANN[i]=c
GNAZ_MARK=[(63,'Helical domain'),(175,'Switch I'),(200,'Switch II'),(228,'Switch III')]
CLS2COL={'sw1':C_SW1,'sw2':C_SW2,'sw3':C_SW3,'hel':C_HELIX,'hd':C_HD}

# RGSZ1 helices (Tesmer, human O76081 numbering): a2-a9 from model DSSP (+259); a1 by homology(approx)
HELICES=[('α1',253,260),('α2',261,267),('α3',270,282),('α4',286,300),('α5',304,319),
         ('α6',331,341),('α7',350,364),('α8',365,370),('α9',372,379)]
RGS_ANN={}
for lab,a,b in HELICES:
    for i in range(a,b+1): RGS_ANN[i]='hel'
RGS_DOMAIN=(262,378)
RGS_MARK=[(a,lab) for lab,a,b in HELICES]

def block(seq, ann, marks, domain=None, per=60):
    """render one protein as monospace lines with a marker row, numbering, colored residues."""
    out=[]
    markpos={p:t for p,t in marks}
    n=len(seq)
    for start in range(0,n,per):
        seg=seq[start:start+per]
        # marker row (clamped so labels fit; suppressed when empty)
        mrow=[' ']*per; has=False
        for p,t in marks:
            if start< p <= start+per:
                col=p-1-start
                if col+len(t)>per: col=per-len(t)
                if col<0: col=0
                for k,ch in enumerate(t):
                    if 0<=col+k<per: mrow[col+k]=ch
                has=True
        if has:
            mrow_html=''.join('&nbsp;' if c==' ' else html.escape(c) for c in mrow)
            out.append('<div class=\"mk\">%s%s</div>'%('&nbsp;'*6, mrow_html))
        # sequence row
        cells=['<span class=\"num\">%4d</span>&nbsp;&nbsp;'%(start+1)]
        for j,aa in enumerate(seg):
            resi=start+j+1
            cls=ann.get(resi)
            style=''
            if cls: style='background:%s;'%CLS2COL[cls]
            if domain and domain[0]<=resi<=domain[1]:
                style+='border-bottom:2px solid #8a8a8a;'
            cells.append('<span style=\"%s\">%s</span>'%(style, html.escape(aa)))
        cells.append('&nbsp;&nbsp;<span class=\"num\">%d</span>'%(min(start+per,n)))
        out.append('<div class=\"seq\">%s</div>'%''.join(cells))
    return '\n'.join(out)

# ---- BLOSUM62 NW for GNAZ vs GNAI1 ----
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
rows=[re.findall(r'-?\d+',r) for r in _bl.strip().split('\n')]
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
    i,j=n,m;A='';C=''
    while i>0 and j>0:
        if F[i][j]==F[i-1][j-1]+sc(a[i-1],b[j-1]):A=a[i-1]+A;C=b[j-1]+C;i-=1;j-=1
        elif F[i][j]==F[i-1][j]+gap:A=a[i-1]+A;C='-'+C;i-=1
        else:A='-'+A;C=b[j-1]+C;j-=1
    while i>0:A=a[i-1]+A;C='-'+C;i-=1
    while j>0:A='-'+A;C=b[j-1]+C;j-=1
    return A,C
alnZ,alnI=nw(GNAZ,GNAI1)

def aln_block(aZ,aI,gnaz_ann,per=60):
    out=[];rz=0;ri=0
    ident=sum(1 for x,y in zip(aZ,aI) if x==y and x!='-')
    for start in range(0,len(aZ),per):
        segZ=aZ[start:start+per];segI=aI[start:start+per]
        z0=rz+1
        # GNAZ row (colored by switch)
        cz=['<span class=\"num\">%4d</span>&nbsp;GNAZ&nbsp;&nbsp;'%z0]
        mid=['&nbsp;'*11]
        ci=[]
        i0=ri+1
        for x,y in zip(segZ,segI):
            zi=rz+1 if x!='-' else rz
            st=''
            if x!='-':
                rz+=1; cls=gnaz_ann.get(rz)
                if cls: st='background:%s;'%CLS2COL[cls]
            cz.append('<span style=\"%s\">%s</span>'%(st, html.escape(x)))
            if y!='-': ri+=1
            mid.append('|' if (x==y and x!='-') else (':' if x!='-' and y!='-' and sc(x,y)>0 else '&nbsp;'))
            ci.append(html.escape(y))
        ci=['<span class=\"num\">%4d</span>&nbsp;GNAI1&nbsp;'%i0]+ci
        out.append('<div class=\"seq\">%s</div>'%''.join(cz))
        out.append('<div class=\"mid\">%s</div>'%''.join(mid))
        out.append('<div class=\"seq\">%s</div>'%''.join(ci))
        out.append('<div class=\"sp\"></div>')
    pct=100.0*ident/sum(1 for c in aZ if c!='-')
    return '\n'.join(out), pct

alnhtml,pct=aln_block(alnZ,alnI,GNAZ_ANN)

leg=lambda col,txt:'<span class=\"chip\" style=\"background:%s\">&nbsp;&nbsp;</span> %s'%(col,txt)
HTML=f'''<div class="wrap">
<h1>Supplementary Figure 1 &mdash; Primary sequences and structural annotations</h1>
<p class="sub">Human sequences from UniProt: GNAZ (P19086, 355 aa), RGS20/RGSZ1 (O76081, 388 aa), GNAI1 (P63096, 354 aa).</p>

<h2>A &nbsp; G&#945;z / GNAZ (P19086) &mdash; switch regions</h2>
<div class="mono">{block(GNAZ,GNAZ_ANN,GNAZ_MARK)}</div>

<h2>B &nbsp; RGSZ1 / RGS20 (O76081) &mdash; RGS domain (underlined, 262&ndash;378) and RGS-box helices (Tesmer &#945;1&ndash;&#945;9)</h2>
<div class="mono">{block(RGS,RGS_ANN,RGS_MARK,domain=RGS_DOMAIN)}</div>

<h2>C &nbsp; GNAZ vs GNAI1 pairwise alignment ({pct:.0f}% identity; &#124; identical, : similar) &mdash; GNAZ switches highlighted</h2>
<div class="mono">{alnhtml}</div>

<div class="legend">
<b>Legend:</b>&nbsp;&nbsp; {leg(C_HD,'G&#945;z helical domain (63&ndash;174)')} &nbsp;&nbsp; {leg(C_SW1,'Switch I (175&ndash;188)')} &nbsp;&nbsp; {leg(C_SW2,'Switch II (200&ndash;220)')} &nbsp;&nbsp; {leg(C_SW3,'Switch III (228&ndash;240)')} &nbsp;&nbsp; {leg(C_HELIX,'RGS-box helix (&#945;1&ndash;&#945;9)')} &nbsp;&nbsp; <span style="border-bottom:2px solid #8a8a8a">underline</span> = annotated RGS domain.
<br><span class="note">Switch ranges (GNAZ) and helix ranges (Tesmer convention, O76081 numbering: &#945;1 253&ndash;260 [approx., N-terminal to the annotated domain], &#945;2 261&ndash;267, &#945;3 270&ndash;282, &#945;4 286&ndash;300, &#945;5 304&ndash;319, &#945;6 331&ndash;341, &#945;7 350&ndash;364, &#945;8 365&ndash;370, &#945;9 372&ndash;379). Helix boundaries from DSSP on the modeled RGS domain mapped to full-length numbering.</span></div>
</div>'''

STYLE='''
*{box-sizing:border-box} body{margin:0;background:#fff}
.wrap{padding:26px 30px;font-family:'Segoe UI',system-ui,sans-serif;color:#111;width:1180px}
h1{font-size:20px;margin:0 0 4px} .sub{font-size:12.5px;color:#555;margin:0 0 14px}
h2{font-size:14px;margin:18px 0 6px;color:#1a1a1a;border-bottom:1px solid #e2e2dc;padding-bottom:3px}
.mono{font-family:'Consolas','Courier New',monospace;font-size:13px;line-height:1.15;white-space:pre}
.mono span{padding:1px 0}
.seq{letter-spacing:1px}
.mk{color:#7a4fd0;font-weight:700;letter-spacing:1px;height:15px}
.mid{color:#888;letter-spacing:1px}
.num{color:#999;font-weight:400;letter-spacing:0}
.sp{height:7px}
.legend{margin-top:16px;font-size:12px;color:#333;border-top:1px solid #e2e2dc;padding-top:10px}
.chip{display:inline-block;width:16px;height:12px;border:1px solid #999;vertical-align:middle}
.note{color:#666;font-size:11px}
'''
open(os.path.join(SCR,'_seqfig.html'),'w',encoding='utf-8').write('<!doctype html><html><head><meta charset="utf-8"><style>%s</style></head><body>%s</body></html>'%(STYLE,HTML))
print('GNAZ %d  RGS20 %d  GNAI1 %d  | GNAZ-GNAI1 identity %.1f%%'%(len(GNAZ),len(RGS),len(GNAI1),pct))
print('wrote _seqfig.html')
