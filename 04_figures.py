"""Figures 1 and 2 (greyscale and colour versions) from results/ files."""
import os, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = 'results/' if os.path.isdir('results') and os.path.exists('results/2_prevalans_sonuclari.csv') else ''
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
OI=['#0072B2','#56B4E9','#009E73','#D55E00','#E69F00','#CC79A7','#000000']
d=pd.read_csv(R+'2_prevalans_sonuclari.csv'); t=d[d.senaryo.str.startswith('S1')]
e=pd.read_csv(R+'5_ekonomik_sonuclar.csv'); s=e[e.senaryo.str.startswith('S1')]
G=[('Erkek 15-18','Men 15–18 y','-','o'),('Erkek 19-64','Men 19–64 y','-','s'),('Erkek 65+','Men ≥65 y','-','^'),
   ('Kadın 15-18','Women 15–18 y','--','o'),('Kadın 19-64','Women 19–64 y','--','s'),('Kadın 65+','Women ≥65 y','--','^'),('Gebe','Pregnant women',':','D')]
for colour in [False, True]:
    fig,ax=plt.subplots(1,2,figsize=(7.2,3.4),sharex=True)
    for k,(th,lab) in enumerate([('lt20','25(OH)D <20 ng/mL (50 nmol/L)'),('lt12','25(OH)D <12 ng/mL (30 nmol/L)')]):
        for i,(g,name,ls,mk) in enumerate(G):
            x=t[t.grup==g].sort_values('doz_mcg_100g')
            c=OI[i] if colour else ('0.45' if g=='Gebe' else 'black')
            ax[k].plot(x.doz_mcg_100g,x[f'prev_{th}_med'],ls=ls,marker=mk,ms=3.5,lw=1.2 if colour else 1,color=c,mfc='white' if 'Kadın' in g else c,label=name)
        ax[k].set_title(('A  ' if k==0 else 'B  ')+lab,loc='left',fontsize=9)
        ax[k].set_xlabel('Fortification level (µg vitamin D$_3$/100 g flour)'); ax[k].set_ylim(0,90 if k==0 else 50); ax[k].set_xlim(0,12.5)
        ax[k].set_xticks([0,1.4,2.8,4,6,8,12]); ax[k].set_xticklabels(['0','1.4','2.8','4','6','8','12'])
        for sp in ['top','right']: ax[k].spines[sp].set_visible(False)
    ax[0].set_ylabel('Prevalence (%)')
    h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=4,frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,0.13,1,1)); n='Figure1'+('_colour' if colour else '')
    fig.savefig(n+'.png',dpi=300); fig.savefig(n+'.tif',dpi=300); plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(7.2,3.3),sharey=True)
    st=[('0','-','o','0'),('0.05','-','s','0.05'),('0.1','--','^','0.10'),('0.25','-.','D','0.25'),('0.5',':','v','0.50')]
    for k,(th,lab) in enumerate([('lt20','Cases of 25(OH)D <20 ng/mL averted'),('lt12','Cases of 25(OH)D <12 ng/mL averted')]):
        for i,(F,ls,mk,name) in enumerate(st):
            ax[k].plot(s.doz,s[f'USD_per_vaka_{th}_F{F}'],ls=ls,marker=mk,ms=3.5,lw=1.2 if colour else 1,color=[OI[0],OI[1],OI[2],OI[4],OI[3]][i] if colour else 'black',mfc='white',label=f'US${name}')
        ax[k].set_title(('A  ' if k==0 else 'B  ')+lab,loc='left',fontsize=9); ax[k].set_ylim(0,13); ax[k].set_xlim(0,12.5)
        ax[k].set_xticks([1.4,2.8,4,6,8,12]); ax[k].set_xticklabels(['1.4','2.8','4','6','8','12'])
        ax[k].set_xlabel('Fortification level (µg vitamin D$_3$/100 g flour)')
        for sp in ['top','right']: ax[k].spines[sp].set_visible(False)
    ax[0].set_ylabel('Cost per case averted (2024 US$)')
    h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=5,frameon=False,fontsize=8,title='Fixed programme cost per capita per year',title_fontsize=8)
    fig.tight_layout(rect=(0,0.15,1,1)); n='Figure2'+('_colour' if colour else '')
    fig.savefig(n+'.png',dpi=300); fig.savefig(n+'.tif',dpi=300); plt.close(fig)
print('figures done')
