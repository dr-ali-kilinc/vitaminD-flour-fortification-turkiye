# Vitamin D fortification of wheat flour in Türkiye: simulation model

Code and results for the manuscript *"Vitamin D fortification of wheat flour in Türkiye: a simulation of effectiveness, safety and costs using national survey data"* (submitted).

## Contents

| File | Description |
|---|---|
| `01_simulation.py` | Fits log-normal baseline 25(OH)D distributions, runs the two-level Monte Carlo simulation for all fortification levels and scenarios, and the worst-case safety analysis |
| `02_safety_by_dose.py` | Proportion with serum 25(OH)D >50 ng/mL by fortification level (main and worst-case assumptions) |
| `03_economic_analysis.py` | Cases averted, fortificant costs and cost per case averted (threshold analysis on fixed programme costs) |
| `results/1_lognormal_uyum.csv` | Observed vs fitted baseline prevalence by population group |
| `results/2_prevalans_sonuclari.csv` | Prevalence of 25(OH)D <20 and <12 ng/mL by group, scenario and fortification level (median, 95% UI) |
| `results/3_guvenlik_sonuclari.csv` | Worst-case safety results (total intake percentiles, % above UL, % above 50/60/100 ng/mL) |
| `results/4_serum_50_uzeri_doz.csv` | % above 50 ng/mL by fortification level |
| `results/5_ekonomik_sonuclar.csv` | Economic results |

Column names and comments in the code are partly in Turkish (e.g. *grup* = group, *senaryo* = scenario, *doz* = fortification level, *onlenen* = averted, *Erkek/Kadın* = men/women, *Gebe* = pregnant, *Emziren* = lactating).

## Data sources

All inputs are published aggregate data and are entered directly in `01_simulation.py`:
- Turkey Nutrition and Health Survey (TBSA) 2017, Ministry of Health, 2019: serum 25(OH)D category distributions, bread and pastry consumption, dietary vitamin D intake, supplement use.
- Turkish Statistical Institute (TurkStat): population 2017, live births 2017.
- Dose–response: restricted cubic spline curves digitised from Dunlop et al., J Nutr 2021;151:2622–35 (Figure 2B and 2D; values in `01_simulation.py`, object `RCS`).
- Costs: WHO/FAO Guidelines on Food Fortification with Micronutrients (2006); wheat flour supply from the Global Fortification Data Exchange.

## How to run

Python 3.12 was used.

```
pip install -r requirements.txt
python 01_simulation.py        # approx. 15–30 min
python 02_safety_by_dose.py
python 03_economic_analysis.py
```

Results are written to `results/`. Random seeds are fixed, so results are reproducible.

## License

MIT License (see `LICENSE`).
