# Vitamin D fortification of wheat flour in Türkiye: simulation model

Code and results for the manuscript *"Vitamin D fortification of wheat flour in Türkiye: a simulation of effectiveness, safety and costs using national survey data"* (submitted).

## Contents

| File | Description |
|---|---|
| `01_simulation.py` | Fits log-normal baseline 25(OH)D distributions, runs the two-level Monte Carlo simulation for all fortification levels and scenarios, and the worst-case safety analysis |
| `02_safety_by_dose.py` | Proportion with serum 25(OH)D >50 ng/mL by fortification level (main and worst-case assumptions) |
| `03_economic_analysis.py` | Cases averted, QALYs, fortificant costs and cost per case averted (threshold analysis on fixed programme costs) |
| `04_figures.py` | Figures 1 and 2 (greyscale and colour) |
| `results/1_lognormal_uyum.csv` | Observed vs fitted baseline prevalence by population group |
| `results/2_prevalans_sonuclari.csv` | Prevalence of 25(OH)D <20 and <12 ng/mL by group, scenario and fortification level (median, 95% UI) |
| `results/3_guvenlik_sonuclari.csv` | Worst-case safety results |
| `results/4_serum_50_uzeri_doz.csv` | % above 50 ng/mL by fortification level |
| `results/5_ekonomik_sonuclar.csv` | Economic results, including QALYs |
| `verification/Dogrulama_Calisma_Kitabi.xlsx` | Independent spreadsheet re-implementation of the main calculations (in Turkish) |

Column names and comments are partly in Turkish (*grup* = group, *senaryo* = scenario, *doz* = fortification level, *onlenen* = averted, *Erkek/Kadın* = men/women, *Gebe* = pregnant, *Emziren* = lactating).

## Data sources

All inputs are published aggregate data entered directly in `01_simulation.py`:
- Turkey Nutrition and Health Survey (TBSA) 2017, Ministry of Health, 2019: serum 25(OH)D categories, design effects for vitamin D categories, bread and pastry consumption, dietary vitamin D intake, supplement use.
- Turkish Statistical Institute (TurkStat): population and live births, 2017.
- Dose–response: restricted cubic spline curves digitised from Dunlop et al., J Nutr 2021;151:2622–35 (object `RCS` in `01_simulation.py`).
- Costs: WHO/FAO Guidelines on Food Fortification with Micronutrients (2006); wheat flour supply from the Global Fortification Data Exchange.

## How to run

Python 3.12 was used.

```
pip install -r requirements.txt
python 01_simulation.py        # approx. 15–30 min
python 02_safety_by_dose.py
python 03_economic_analysis.py
python 04_figures.py
```

Random seeds are fixed, so results are reproducible.

## Version history

- **v1.1.0**: Sampling uncertainty now uses the design effect reported in TBSA 2017 for serum vitamin D categories (DEFF 3.07) instead of an estimate derived from confidence intervals; 19–64-year sub-bands pooled with 2017 population weights instead of sample sizes; QALY outputs (including a scenario restricted to cases <12 ng/mL) added to the economic results; figure script and verification workbook added.
- **v1.0.0**: Initial version.

## License

MIT License (see `LICENSE`).
