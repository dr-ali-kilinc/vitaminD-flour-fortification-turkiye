"""
Ekonomik analiz: önlenen vaka sayısı ve önlenen vaka başına maliyet (eşik yaklaşımı)
Maliyet = sabit program maliyeti (F, kişi başı/yıl; eşik analizi) + doza bağlı D3 hammadde maliyeti
Hammadde: WHO/FAO Guidelines on Food Fortification (Allen et al. 2006), "Estimated costs of selected fortificants":
  suda çözünür kuru D3, içerik %0.25, 33 USD/kg (2006 USD), %20 overage
Un hacmi: GFDx Türkiye buğday unu arzı 337 g/kişi/gün (tüm nüfus) -> tüm buğday unu zenginleştirilir varsayımı (üst sınır)
Nüfus: TÜİK ADNKS 2017; gebe = 2017 canlı doğum (1.300.258) x 0.75
"""
import numpy as np, pandas as pd
CPI = 1.5553               # BLS CPI Inflation Calculator, Ocak 2006 -> Ocak 2024
PRICE, CONTENT, OVERAGE = 33*CPI, 0.0025, 1.20
POP_TOTAL, FLOUR_G = 80_810_525, 337
births = 1_300_258; preg = births*0.75
POP = {'Erkek 15-18':2_680_834,'Erkek 19-64':25_051_767,'Erkek 65+':3_033_433,
       'Kadın 15-18':2_540_435,'Kadın 19-64':24_608_616-preg,'Kadın 65+':3_861_952,'Gebe':preg}

def ingredient_cost(dose):   # USD/yıl, tüm ülke
    flour_t = FLOUR_G*1e-6*365*POP_TOTAL
    d3_kg = dose*10*flour_t/1e6               # dose mcg/100g = 10*dose mg/ton
    return d3_kg/CONTENT*OVERAGE*PRICE

d = pd.read_csv('results/2_prevalans_sonuclari.csv')
rows=[]
for sc in ['S1','S9','S6','S7']:
    for dose in sorted(d.doz_mcg_100g.unique()):
        if dose==0: continue
        sub = d[(d.senaryo.str.startswith(sc))&(d.doz_mcg_100g==dose)]
        if sub.empty: continue
        r={'senaryo':sub.senaryo.iloc[0],'doz':dose}
        for th in ['lt20','lt12']:
            # delta negatif: alt95 (en negatif) -> en çok önlenen vaka; etiketler buna göre çevrildi
            for q_out,q_in in [('med','med'),('alt95','ust95'),('ust95','alt95')]:
                r[f'onlenen_{th}_{q_out}'] = -sum(sub[sub.grup==g][f'delta_{th}_{q_in}'].iloc[0]/100*POP[g] for g in POP)
        ic = ingredient_cost(dose)
        r['hammadde_USD_yil']=ic; r['hammadde_USD_kisi_yil']=ic/POP_TOTAL
        for F in [0, 0.05, 0.10, 0.25, 0.50]:
            tot = ic + F*POP_TOTAL
            r[f'USD_per_vaka_lt20_F{F}'] = tot/r['onlenen_lt20_med']
            r[f'USD_per_vaka_lt12_F{F}'] = tot/r['onlenen_lt12_med']
        rows.append(r)
out=pd.DataFrame(rows); out.to_csv('results/5_ekonomik_sonuclar.csv',index=False,encoding='utf-8-sig')
pd.set_option('display.width',250); pd.set_option('display.float_format',lambda x:f'{x:,.2f}')
s=out[out.senaryo.str.startswith('S1')]
print(s[['doz','onlenen_lt20_med','onlenen_lt20_alt95','onlenen_lt20_ust95','onlenen_lt12_med','hammadde_USD_yil','hammadde_USD_kisi_yil']].to_string(index=False))
print(s[['doz']+[c for c in out.columns if 'per_vaka' in c]].to_string(index=False))
print(out[out.senaryo.str[:2].isin(['S9','S6','S7'])&out.doz.isin([2.8,4,8])][['senaryo','doz','onlenen_lt20_med','onlenen_lt12_med','USD_per_vaka_lt20_F0.1','USD_per_vaka_lt12_F0.1']].to_string(index=False))
