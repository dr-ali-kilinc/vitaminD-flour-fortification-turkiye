"""
Proportion of individuals with serum 25(OH)D >50 ng/mL (125 nmol/L) by fortification level,
under main-analysis assumptions (Dunlop adult curve, flour content 0.74, retention 0.85) and
worst-case assumptions (bread + pastries, upper flour content, retention 0.90, 0.85 ng/mL per ug/day).
Requires 01_simulation.py in the same folder (re-uses its input data and fitted distributions).
"""
import numpy as np, pandas as pd, importlib.util, sys, io, contextlib
src=open('01_simulation.py').read()
src=src.split('# ---------------- 4. Senaryolar')[0].replace('N_OUTER, N_INNER = 20251014, 1000, 20000','N_OUTER, N_INNER = 20251014, 2, 20000')
ns={}; exec(src, ns)
fits={f['grup']:f for f in ns['fits']}; INT=ns['INTAKE']; ln=ns['ln_draw']; rcs=ns['rcs_increment']
rs=np.random.default_rng(7); N=200000; rows=[]
for g,f in fits.items():
    b_m,b_sd,c_m,c_sd,*_=INT[g]
    base=np.exp(f['mu_log']+f['sigma_log']*rs.standard_normal(N))
    bread=ln(b_m,b_sd,rs.standard_normal(N)); corek=ln(c_m,c_sd,rs.standard_normal(N))
    r={'grup':g,'bazal':round(np.mean(base>50)*100,2)}
    for d in [2,4,6,8,12]:
        main=d/100*bread*0.74*0.85
        worst=d/100*(bread*0.79+corek*0.65)*0.90
        r[f'ana_{d}']=round(np.mean(base+rcs('yetiskin',main,0)>50)*100,2)
        r[f'enkotu_{d}']=round(np.mean(base+0.85*worst>50)*100,2)
        r[f'ortalama_alim_{d}']=round(main.mean(),1)
    rows.append(r)
df=pd.DataFrame(rows); df.to_csv('results/4_serum_50_uzeri_doz.csv',index=False,encoding='utf-8-sig')
pd.set_option('display.width',250); print(df.to_string(index=False))
