#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
GA=os.path.dirname(os.path.abspath(__file__)); OUT=r"C:\Users\richi\Documents\RGSZ1-GNAZ Project"
L=json.load(open(os.path.join(GA,'labels.json')))
M=np.load(os.path.join(GA,'gal_idmat.npy')); R=np.load(os.path.join(GA,'rgs_idmat.npy'))
plt.rcParams.update({'font.family':'DejaVu Sans'})
cmap=LinearSegmentedColormap.from_list('id',['#f7fbff','#c6dbef','#6baed6','#2171b5','#08306b'])

# ---- Galpha 5x5 ----
gn=L['gnames']; n=len(gn)
fig,ax=plt.subplots(figsize=(5.6,5.0))
im=ax.imshow(M,cmap=cmap,vmin=55,vmax=100)
for i in range(n):
    for j in range(n):
        ax.text(j,i,'%.0f'%M[i,j] if i!=j else '—',ha='center',va='center',
                fontsize=12,color='white' if M[i,j]>82 else '#222',fontweight='bold' if i==j else 'normal')
ax.set_xticks(range(n)); ax.set_xticklabels(gn,fontsize=11)
ax.set_yticks(range(n)); ax.set_yticklabels(gn,fontsize=11)
ax.set_title('Human Gα subunits — pairwise % identity',fontsize=13,fontweight='bold',pad=10)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.set_label('% identity',fontsize=10); cb.outline.set_visible(False)
fig.text(0.01,0.01,'Global Needleman–Wunsch, BLOSUM62; full-length UniProt sequences. Identity over aligned (ungapped) columns.',fontsize=7.5,color='#666')
fig.tight_layout(rect=[0,0.03,1,1]); fig.savefig(OUT+r"\galpha_identity_heatmap.png",dpi=200,bbox_inches='tight'); plt.close(fig)

# ---- RGS 17x17 ----
rn=L['rn']; fam=L['fam']; K=len(rn)
lab=['%s'%rn[i] for i in range(K)]
fig,ax=plt.subplots(figsize=(10.4,9.4))
im=ax.imshow(R,cmap=cmap,vmin=28,vmax=100)
for i in range(K):
    for j in range(K):
        v=R[i,j]
        ax.text(j,i,'%.0f'%v if i!=j else '—',ha='center',va='center',
                fontsize=8.2,color='white' if v>72 else '#333')
# family boundaries
bnds=[k for k in range(1,K) if fam[k]!=fam[k-1]]
for b in bnds:
    ax.axhline(b-0.5,color='#b22',lw=1.6); ax.axvline(b-0.5,color='#b22',lw=1.6)
ax.set_xticks(range(K)); ax.set_xticklabels(lab,fontsize=9,rotation=90)
ax.set_yticks(range(K)); ax.set_yticklabels(lab,fontsize=9)
# family band labels
segs=[0]+bnds+[K]; fams=[]
for a,b in zip(segs[:-1],segs[1:]): fams.append((fam[a],(a+b-1)/2))
for fa,mid in fams:
    ax.text(mid,-0.85,fa,ha='center',va='bottom',fontsize=12,fontweight='bold',color='#b22')
    ax.text(-1.7,mid,fa,ha='center',va='center',fontsize=12,fontweight='bold',color='#b22',rotation=90)
ax.set_title('Human RGS domains (R4 / R7 / RZ) — pairwise % identity',fontsize=14,fontweight='bold',pad=44)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.02); cb.set_label('% identity',fontsize=10); cb.outline.set_visible(False)
fig.text(0.01,0.005,'RGS-box residues only (UniProt domain boundaries); global Needleman–Wunsch, BLOSUM62; identity over aligned (ungapped) columns. Red lines separate subfamilies.',fontsize=8,color='#666')
fig.tight_layout(rect=[0,0.02,1,1]); fig.savefig(OUT+r"\rgs_identity_heatmap.png",dpi=200,bbox_inches='tight'); plt.close(fig)
print('wrote galpha_identity_heatmap.png and rgs_identity_heatmap.png')
