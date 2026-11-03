# FOI June 2013 TNT and TNT/Al shots: 1D free-air stage

Validation data for the products and Al afterburning (`report/h0_dataset_survey.md`,
candidate 9): E. Fedina, "TNT/Aluminium afterburning in air blasts",
FOI-R--3913--SE (2014), Chapter 4, and Fedina, PhD thesis, Lund University (2017).

## Experiment

- 1 kg spheres of TNT and TNT/Al 80/20 (A100 Al, 42 um), 2 kg TNT; 30 g (45 g
  for 2 kg) m/46 booster in a cavity at the top. Hung at HOB 0.15, 0.5, 1.0 and
  2.0 m, one TNT and one TNT/Al shot each (TNT at 0.5 m twice).
- Ground-flush PCB gauges P1-P3 at 2, 4, 6 m from the point below the charge
  (1 MHz); pyrometer aimed about 1.5 m up at 2 m radius (650-1800 C).
- `foi2013_data.py`: arrival, peak, positive impulse and duration at P1-P3 and
  the pyrometer maximum of each shot, from the raw files with
  `extract_shots.py` (the raw data stay outside the repository; set
  `FOI2013_DATA` to the `Measurements/TXT` folder). Shot N is the file
  `Katja SkottN`; the TNT peaks agree with Table 4.2 of the report.
- Repeatability (shots 2 and 9): 3% in P2/P3 impulse, 10-11% in peak. Al effect
  at 4 and 6 m: +15 to +27% in positive impulse at every HOB.

## Case

The gauges are on the ground, so the comparison needs the 2D axisymmetric
case with a ground plane. This 1D spherical case is its free-air stage: it
runs until (and beyond) the shock reaches the ground, and its fields are mapped
into 2D.

```
python3 setup_foi2013.py --shot 4 --out <dir> [--closure laminar|pasrd]
cd <dir> && ./Allrun
```

- Built with `build_case()` of `../Pontalier_NM_1D_spherical/setup_pontalier.py`
  (blastEulerFoam, products -> gas transfer below 100 bar, 28-species gas with
  equilibrium above 2500 K, Phase D Al models).
- Products: EuDetCode products-only fit (BKW) of the main charge with the booster
  lumped in (`utilities/thermo_fitter/foi2013_charges.py`,
  `eudet_foi2013/shot<N>`); shot 4: D = 6697 m/s, p_CJ = 16.9 GPa, 17.5%
  graphite in the transferred products. TNT/Al: 90% of the Al as the Al phase.
- Assumed: TNT 1630 kg/m3, booster PETN/oil 86/14 at 1500 kg/m3, ambient
  101325 Pa and 293.15 K (not recorded).
- Domain 4 m (0.5 mm to 0.1 m, graded to 15 mm), probes at the slant distances
  of the gauges inside the domain.
