"""
Buğday unu D3 vitamini fortifikasyonu - birey düzeyinde iki katmanlı Monte Carlo simülasyonu
Veri kaynağı: TBSA 2017 (yayın 2019) raporu
  - Serum 25(OH)D dağılım dilimleri: Tablo 3.73 (15+ yaş), Tablo 3.83 (gebe)
  - Ekmek ve "çörek vb." tüketimi (g/gün, ort ± SS): Tablo 2.55, 2.56, 2.60
  - Diyetle D vitamini alımı (mcg/gün, ort ± SS): Tablo 2.53, 2.54, 2.59
  - D vitamini takviyesi kullanım oranı: Tablo 2.48
Dış katman: parametre belirsizliği (N_OUTER iterasyon)
İç katman: popülasyon değişkenliği (N_INNER sanal birey / grup)
"""
import numpy as np, pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize

SEED, N_OUTER, N_INNER = 20251014, 1000, 20000
rng = np.random.default_rng(SEED)
DOSES = [0, 1.4, 2, 2.8, 4, 6, 8, 12]   # 1.4: GCC GSO/FDS 194/2014; 2.8: 2x GCC (Cashman & O'Dea 2019)                       # mcg D3 / 100 g un
EDGES = np.array([0, 10, 20, 30, 80, np.inf])

# ---------------- 1. TBSA girdileri ----------------
# (N, ortalama, %GA alt, %GA üst, [<10, 10-19, 20-29, 30-79, >=80] ağırlıklı %)
SERUM = {
 'Erkek 15-18': [(207,23.5,22.0,24.9,[2.9,35.7,37.1,24.3,0])],
 'Erkek 19-64': [(752,20.4,19.7,21.2,[7.5,43.4,35.2,13.7,0.2]),(2089,19.9,19.5,20.4,[7.1,47.3,32.8,12.5,0.2]),(1102,21.2,20.5,21.9,[6.8,41.4,37.3,14.1,0.3])],
 'Erkek 65+':   [(847,21.4,20.5,22.3,[8.6,42.2,31.5,17.4,0.2])],
 'Kadın 15-18': [(237,13.4,12.5,14.4,[30.5,53.5,14.0,2.0,0])],
 'Kadın 19-64': [(990,13.8,13.2,14.5,[36.3,44.9,15.4,3.2,0.1]),(2346,15.9,15.3,16.5,[27.3,47.2,18.0,6.9,0.5]),(1378,17.3,16.5,18.1,[23.3,47.4,19.1,9.9,0.4])],
 'Kadın 65+':   [(1119,18.2,16.9,19.5,[31.5,34.7,18.7,14.6,0.5])],
 'Gebe':        [(137,14.9,13.3,16.6,[24.3,47.2,23.6,4.9,0])],
}
# (ekmek ort, SS, çörek ort, SS, diyet D ort, SS, takviye oranı, takviye dozu mcg)
INTAKE = {
 'Erkek 15-18': (254.8,159.39, 35.3,67.87, 2.4,3.72, 0.004,25),
 'Erkek 19-64': (227.8,143.72, 19.8,41.93, 3.7,8.77, 0.022,25),
 'Erkek 65+':   (188.5,125.35,  9.6,26.44, 3.5,9.78, 0.028,25),
 'Kadın 15-18': (138.4, 91.10, 36.8,46.34, 2.0,4.92, 0.004,25),
 'Kadın 19-64': (134.5, 96.70, 17.0,32.63, 3.1,13.02,0.022,25),
 'Kadın 65+':   (127.2, 86.04,  9.0,29.81, 2.7,7.73, 0.028,25),
 'Gebe':        (151.3, 84.17, 22.1,37.83, 4.3,21.97,1.000,30),  # program: tüm gebeler 1200 IU (Genelge 2011/34; en kötü durum)
 'Emziren':     (172.2,136.06, 20.9,34.94, 2.8,5.33, 1.000,30),  # serum verisi yok; program: doğum sonrası 6 ay 1200 IU (Genelge 2011/34), tümüne atandı (en kötü durum)
}

# ---------------- 2. Log-normal'in dilimlere uydurulması ----------------
def fit_lognormal(p, n):
    def nll(th):
        mu, s = th[0], np.exp(th[1])
        c = norm.cdf((np.log(np.maximum(EDGES,1e-9))-mu)/s); c[0], c[-1] = 0, 1
        return -np.sum(n*p*np.log(np.clip(np.diff(c),1e-12,1)))
    o = minimize(nll, [np.log(15), -0.6], method='Nelder-Mead', options={'xatol':1e-5,'fatol':1e-5})
    return o.x[0], np.exp(o.x[1])

# 19-64 yaş alt bantları (19-30, 31-50, 51-64) TÜİK 2017 nüfusuyla ağırlıklandırılır
POPW = {'Erkek 19-64': [7766249, 11777664, 5507854], 'Kadın 19-64': [7486506, 11551635, 5570475]}
# Desen etkisi: TBSA 2017 Ek, "Laboratuvar Değişkenleri İçin Hesaplanan Desen Etkileri" (s. 437),
# D vitamini dilimleri DEFT = 1.722, 1.804, 1.730, 1.633, 1.231 -> DEFF = DEFT^2, dilim yüzdeleriyle ağırlıklı ortalama
DEFT = np.array([1.722, 1.804, 1.730, 1.633, 1.231]); DEFT_P = np.array([18.2, 44.7, 25.9, 10.9, 0.3])
DEFF = float(np.sum(DEFT**2*DEFT_P)/DEFT_P.sum())   # = 3.065

def pooled(g, rows):
    N = sum(r[0] for r in rows)                       # toplam örneklem (etkin N için)
    w = np.array(POPW.get(g, [r[0] for r in rows]), float)
    p = sum(wi*np.array(r[4],float) for wi, r in zip(w, rows))/w.sum()
    return N, p/p.sum()

fits, param_draws = [], {}
for g, rows in SERUM.items():
    N, p = pooled(g, rows)
    mu, s = fit_lognormal(p, N)
    sd_fit = np.exp(mu+s*s/2)*np.sqrt(np.exp(s*s)-1)
    deff = DEFF; n_eff = N/deff
    # Örnekleme belirsizliği: etkin N ile Dirichlet yeniden örnekleme -> yeniden uydurma
    d = rng.dirichlet(p*n_eff + 0.5, size=N_OUTER)
    param_draws[g] = np.array([fit_lognormal(x/x.sum(), n_eff) for x in d])
    c = norm.cdf((np.log(np.maximum(EDGES,1e-9))-mu)/s); c[0], c[-1] = 0, 1
    fp = np.diff(c)*100
    fits.append(dict(grup=g, N=N, tasarim_etkisi=round(deff,3), mu_log=round(mu,4), sigma_log=round(s,4),
        medyan=round(np.exp(mu),1), ortalama_fit=round(np.exp(mu+s*s/2),1), SS_fit=round(sd_fit,1),
        gozlenen_lt10=round(p[0]*100,1), fit_lt10=round(fp[0],1),
        gozlenen_lt20=round((p[0]+p[1])*100,1), fit_lt20=round(fp[0]+fp[1],1),
        fit_lt12=round(norm.cdf((np.log(12)-mu)/s)*100,1)))
pd.DataFrame(fits).to_csv('results/1_lognormal_uyum.csv', index=False, encoding='utf-8-sig')

def ln_draw(mean, sd, z):
    s = np.sqrt(np.log(1+(sd/mean)**2)); return np.exp(np.log(mean)-s*s/2 + s*z)

# ---------------- 3. Doz-yanıt fonksiyonları ----------------
# Dunlop E, Kiely ME, James AP, Singh T, Pham NM, Black LJ. J Nutr 2021;151:2622-35. doi:10.1093/jn/nxab180
# Şekil 2'den sayısallaştırma (RCS, düğümler 0/5/20 mcg/gün); y = 25(OH)D ortalama fark, nmol/L
DR_X = np.array([0, 5, 10, 15, 20, 25, 30, 40])
RCS = {  # (orta, alt %95, üst %95)
 'yetiskin': (np.array([0,11.0,19.0,23.5,26.0,27.5,29.0,32.5]),
              np.array([0, 8.0,14.0,18.5,21.0,22.2,22.3,22.3]),
              np.array([0,13.5,22.5,27.0,30.0,33.0,36.5,43.5])),
 'D3':       (np.array([0,15.0,24.0,28.5,30.0,30.0,30.0,30.0]),   # 20 mcg sonrası düşüş biyolojik olarak
              np.array([0,10.0,18.0,22.5,24.0,24.0,24.0,24.0]),   # anlamlı değil -> plato kabul edildi
              np.array([0,20.0,30.0,35.0,36.0,37.0,38.5,44.0])),
}
NMOL_PER_NG = 2.496
def rcs_increment(name, intake, z):
    mid, lo, hi = RCS[name]
    y = mid + (z*(hi-mid)/1.96 if z >= 0 else z*(mid-lo)/1.96)
    y = np.maximum.accumulate(np.maximum(y, 0))            # monoton artan
    tail = (y[-1]-y[-2])/(DR_X[-1]-DR_X[-2])                # 40 mcg üstü: son eğimle doğrusal
    inc = np.interp(intake, DR_X, y) + np.where(intake > 40, (intake-40)*tail, 0)
    return inc / NMOL_PER_NG                                # ng/mL

# ---------------- 4. Senaryolar ----------------
SCEN = {
 'S1 Ana: ekmek, %100, Dunlop yetişkin eğrisi':      dict(corek=False, cov=1.00, dr=('rcs','yetiskin')),
 'S2 Ekmek + çörek, %100, yetişkin eğrisi':          dict(corek=True,  cov=1.00, dr=('rcs','yetiskin')),
 'S3 Ekmek, %90, yetişkin eğrisi':                   dict(corek=False, cov=0.90, dr=('rcs','yetiskin')),
 'S4 Ekmek, %75, yetişkin eğrisi':                   dict(corek=False, cov=0.75, dr=('rcs','yetiskin')),
 'S5 Ekmek, %100, Dunlop D3 eğrisi':                 dict(corek=False, cov=1.00, dr=('rcs','D3')),
 'S6 Doğrusal 2 nmol/L/100 IU (alt sınır)':          dict(corek=False, cov=1.00, dr=('lin',0.32,0.04)),
 'S7 Doğrusal 1.2 nmol/L/mcg (Black 2012)':          dict(corek=False, cov=1.00, dr=('lin',0.48,0.06)),
 'S8 Doğrusal 1.7 nmol/L/mcg (Cashman 2021 IPD)':    dict(corek=False, cov=1.00, dr=('lin',0.68,0.08)),
 # Nikooyeh & Neyestani BJN 2022: bazal <50 nmol/L alt grubunda 3 nmol/L/100 IU = 0.48 ng/mL/mcg
 'S9 Fırında yüksek kayıp (muhafaza %50-70)':         dict(corek=False, cov=1.00, dr=('rcs','yetiskin'), ret=(0.50,0.70)),
}
# Belirsiz parametreler (dış katman); doz-yanıt belirsizliği tek z ile (eğri GA'sı / doğrusal eğim SS)
# un oranı ekmek: U(0.69,0.79) - kuru madde dengesi: TÜRKOMP beyaz ekmek (nem 31.64, tuz 0.858 g/100 g) ve TGK Ekmek Tebliği (nem en çok %38); çörek: U(0.45,0.65); fırında muhafaza: U(0.80,0.90)
# doz-yanıt eğimi: ng/mL artış / 1 mcg/gün ek alım ~ Normal (0.32 = 0.8 ng/mL / 100 IU)

out = []
groups = [g for g in SERUM]
for gi, g in enumerate(groups):
    b_m,b_sd,c_m,c_sd,*_ = INTAKE[g]
    rs = np.random.default_rng(SEED + 100*(gi+1))
    res = {(sc,d): {'lt20':[], 'lt12':[], 'mean':[]} for sc in SCEN for d in DOSES}
    for k in range(N_OUTER):
        mu, s = param_draws[g][k]
        zb, z1, z2, u = rs.standard_normal(N_INNER), rs.standard_normal(N_INNER), rs.standard_normal(N_INNER), rs.random(N_INNER)
        base = np.exp(mu + s*zb)
        bread, corek = ln_draw(b_m,b_sd,z1), ln_draw(c_m,c_sd,z2)
        ff_b, ff_c, u_ret = rs.uniform(.69,.79), rs.uniform(.45,.65), rs.random()
        zs = rs.standard_normal()
        for sc, P in SCEN.items():
            lo_r, hi_r = P.get('ret', (0.80, 0.90))
            ret = lo_r + (hi_r-lo_r)*u_ret
            flour = bread*ff_b + (corek*ff_c if P['corek'] else 0)
            covered = (u < P['cov'])
            for d in DOSES:
                intake = d/100*flour*ret*covered            # mcg/gün ek D3
                if P['dr'][0] == 'rcs':
                    serum = base + rcs_increment(P['dr'][1], intake, zs)
                else:
                    serum = base + max(P['dr'][1] + P['dr'][2]*zs, 0)*intake
                r = res[(sc,d)]
                r['lt20'].append(np.mean(serum<20)*100); r['lt12'].append(np.mean(serum<12)*100); r['mean'].append(serum.mean())
    for (sc,d), r in res.items():
        b20 = np.array(res[(sc,0)]['lt20']); b12 = np.array(res[(sc,0)]['lt12'])
        a20, a12, m = map(np.array, (r['lt20'], r['lt12'], r['mean']))
        q = lambda x: np.percentile(x,[50,2.5,97.5])
        out.append(dict(grup=g, senaryo=sc, doz_mcg_100g=d,
            **{f'{k}_{t}':round(v,2) for k,x in [('serum_ort',m),('prev_lt20',a20),('prev_lt12',a12),
                ('delta_lt20',a20-b20),('delta_lt12',a12-b12)] for t,v in zip(['med','alt95','ust95'],q(x))}))
pd.DataFrame(out).to_csv('results/2_prevalans_sonuclari.csv', index=False, encoding='utf-8-sig')

# ---------------- 4. Güvenlik (en kötü durum) ----------------
# 12 mcg/100 g, ekmek + çörek, %100 kapsama, üst sınır un oranı ve muhafaza, takviye dahil
saf = []
rs = np.random.default_rng(SEED+1)
N_S = 200000
for g,(b_m,b_sd,c_m,c_sd,dm,dsd,sup_p,sup_d) in INTAKE.items():
    bread, corek, diet = ln_draw(b_m,b_sd,rs.standard_normal(N_S)), ln_draw(c_m,c_sd,rs.standard_normal(N_S)), ln_draw(dm,dsd,rs.standard_normal(N_S))
    fort = 12/100*(bread*0.79 + corek*0.65)*0.90
    sup = (rs.random(N_S) < sup_p)*sup_d
    total = diet + fort + sup
    row = dict(grup=g, ekmek_P97_5_g=round(np.percentile(bread,97.5)),
        fortifikasyon_P50=round(np.median(fort),1), fortifikasyon_P97_5=round(np.percentile(fort,97.5),1),
        toplam_P97_5=round(np.percentile(total,97.5),1), toplam_P99=round(np.percentile(total,99),1),
        UL_ustu_yuzde=round(np.mean(total>100)*100,3),
        UL_ustu_yuzde_fortifikasyonsuz=round(np.mean((diet+sup)>100)*100,3))
    if g in SERUM:
        mu,s = [f for f in fits if f['grup']==g][0]['mu_log'], [f for f in fits if f['grup']==g][0]['sigma_log']
        base = np.exp(mu + s*rs.standard_normal(N_S))
        serum = base + 0.85*fort                               # en yüksek eğim (S6): muhafazakâr
        # Eşikler: >50 ng/mL (125 nmol/L, IOM "olası zarar" endişesi), >100 ng/mL (250 nmol/L, klasik toksisite)
        # Kaynak: McKenna MJ, Murray BF. Endocr Connect 2013;2:87-95 (IOM 2011'e atıfla)
        row.update(serum_gt50_yuzde=round(np.mean(serum>50)*100,2), serum_gt50_yuzde_fortifikasyonsuz=round(np.mean(base>50)*100,2))
        # IOM endişe aralığının üst ucu 150 nmol/L = 60 ng/mL (Cashman et al. Eur J Nutr 2021;60:939-59)
        row.update(serum_gt60_yuzde=round(np.mean(serum>60)*100,2), serum_gt60_yuzde_fortifikasyonsuz=round(np.mean(base>60)*100,2))
        row.update(serum_P99_9=round(np.percentile(serum,99.9),1), serum_gt100_yuzde=round(np.mean(serum>100)*100,3),
                   serum_P99_9_fortifikasyonsuz=round(np.percentile(base,99.9),1),
                   serum_gt100_yuzde_fortifikasyonsuz=round(np.mean(base>100)*100,3))
    saf.append(row)
pd.DataFrame(saf).to_csv('results/3_guvenlik_sonuclari.csv', index=False, encoding='utf-8-sig')
print('Tamamlandı')
