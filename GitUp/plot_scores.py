#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
OUT=r"C:\Users\richi\Documents\RGSZ1-GNAZ Project"
d=json.load(open(r"C:\Users\richi\Downloads\scores.json"))
MODELS=['1','2','3']
CH=[(0,355,'GNAZ'),(355,472,'RGSZ1'),(472,504,'GTP'),(504,505,'Mg')]  # 0-indexed token ranges
BND=[355,472,504]           # divider positions
CHNAMES=['GNAZ','RGSZ1','GTP','Mg']
MCOL=['#2166ac','#d6604d','#1a9850']  # model 1,2,3
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})

# ---------- 1) pLDDT ----------
fig,ax=plt.subplots(figsize=(11,4.2))
# confidence bands
for lo,hi,c in [(90,100,'#d7ecd9'),(70,90,'#eef7d9'),(50,70,'#fdecc8'),(0,50,'#fbdcd4')]:
    ax.axhspan(lo,hi,color=c,alpha=.55,zorder=0,lw=0)
for m,col in zip(MODELS,MCOL):
    pl=np.array(d[m]['plddt'])*100
    ax.plot(np.arange(1,len(pl)+1),pl,color=col,lw=1.1,label='Model %s (complex pLDDT %.1f)'%(m,d[m]['complex_plddt']*100))
for b in BND: ax.axvline(b+0.5,color='#444',lw=.8,ls='--',alpha=.7)
for a,b,nm in CH: ax.text((a+b)/2+0.5,1.02,nm,ha='center',va='bottom',fontsize=9,fontweight='bold',color='#333',transform=ax.get_xaxis_transform())
ax.set_xlim(1,505); ax.set_ylim(0,100); ax.set_xlabel('Token (residue for chains; per-atom for GTP/Mg)')
ax.set_ylabel('pLDDT'); ax.set_title('Per-token pLDDT (AlphaFold3 / Boltz-2, 3 models)',fontweight='bold',pad=26)
ax.legend(loc='lower right',fontsize=8,framealpha=.9)
ax.margins(x=0); fig.tight_layout(); fig.savefig(OUT+r"\pLDDT_plot.png",dpi=150,bbox_inches='tight'); plt.close(fig)

# ---------- helper for square matrix panels ----------
def matrix_fig(key,cmap,vmax,label,fname,vmin=0):
    fig,axs=plt.subplots(1,3,figsize=(15,5.4))
    for k,(m,ax) in enumerate(zip(MODELS,axs)):
        M=np.array(d[m][key])
        im=ax.imshow(M,cmap=cmap,vmin=vmin,vmax=vmax,origin='upper',interpolation='nearest')
        for b in BND:
            ax.axvline(b-0.5,color='k',lw=.6); ax.axhline(b-0.5,color='k',lw=.6)
        mids=[(a+b)/2 for a,b,_ in CH]
        ax.set_xticks(mids); ax.set_xticklabels(CHNAMES,fontsize=8)
        ax.set_yticks(mids); ax.set_yticklabels(CHNAMES,fontsize=8,rotation=90,va='center')
        ax.set_title('Model %s'%m,fontsize=11,fontweight='bold')
        ax.set_xlabel('Aligned residue');
        if k==0: ax.set_ylabel('Scored residue')
    cb=fig.colorbar(im,ax=axs,fraction=0.025,pad=0.02); cb.set_label(label)
    fig.suptitle(label+'  (505 tokens: GNAZ | RGSZ1 | GTP | Mg)',fontweight='bold',y=0.99)
    fig.savefig(OUT+'\\'+fname,dpi=150,bbox_inches='tight'); plt.close(fig)

# ---------- 2) PAE ----------
matrix_fig('pae','Greens_r',31.75,'Predicted aligned error (Å)','PAE_plot.png')
# ---------- 3) PDE ----------
pdemax=max(np.array(d[m]['pde']).max() for m in MODELS)
matrix_fig('pde','Blues_r',float(np.ceil(pdemax)),'Predicted distance error (Å)','PDE_plot.png')

# ---------- 4) per-chain (pairwise) ipTM ----------
fig,axs=plt.subplots(1,3,figsize=(13,4.6))
for m,ax in zip(MODELS,axs):
    P=d[m]['pair_chains_iptm']
    M=np.array([[P[str(i)][str(j)] for j in range(4)] for i in range(4)])
    im=ax.imshow(M,cmap='viridis',vmin=0,vmax=1)
    for i in range(4):
        for j in range(4):
            ax.text(j,i,'%.2f'%M[i,j],ha='center',va='center',
                    color='white' if M[i,j]<0.6 else 'black',fontsize=9)
    ax.set_xticks(range(4)); ax.set_xticklabels(CHNAMES,fontsize=8)
    ax.set_yticks(range(4)); ax.set_yticklabels(CHNAMES,fontsize=8)
    ax.set_title('Model %s (ipTM %.3f)'%(m,d[m]['iptm']),fontsize=10,fontweight='bold')
    ax.set_xlabel('Chain j');
    if m=='1': ax.set_ylabel('Chain i')
cb=fig.colorbar(im,ax=axs,fraction=0.022,pad=0.02); cb.set_label('pairwise ipTM')
fig.suptitle('Per-chain (pairwise) ipTM',fontweight='bold',y=1.02)
fig.savefig(OUT+r"\ipTM_plot.png",dpi=150,bbox_inches='tight'); plt.close(fig)
print('wrote pLDDT_plot.png, PAE_plot.png, PDE_plot.png, ipTM_plot.png to project folder')
print('PDE max = %.2f'%pdemax)
