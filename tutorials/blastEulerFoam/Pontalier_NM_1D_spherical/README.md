# 1D spherical gelled-nitromethane charges of Pontalier et al. (2021)

Validation case H1 (`report/h0_dataset_survey.md`): free-field blast from
about 1 kg spheres of gelled nitromethane with embedded Al or glass
particles, Q. Pontalier, J. Loiseau, S. Goroshin, F. Zhang, D. L. Frost,
"Blast enhancement from metalized explosives", Shock Waves 31 (2021)
203-230, doi:10.1007/s00193-021-00994-z.

## Experiment

- Gelled NM (95.6 NM / 4 PMMA / 0.4 glass microballoons) in a 12.7 cm glass
  bulb, with Al (13, 54, 114 um; 10-50% by mass), Mg or glass beads (66 um);
  central C-4 booster of 86-93 g. Height of burst 1.52 m, altitude 700 m
  (p0 about 93 kPa).
- Gauges at 1.0, 1.5 and 2.0 m; the ground reflection reaches the 2.0 m gauge
  after the positive phase. One shot per configuration.
- `pontalier_data.py`: arrival times (Table 3) and the 2.0 m peak
  overpressure from the gauge (Table 5) and from the video (Table 7).

## Case

`setup_pontalier.py` builds the case of one trial from
`../ChargeCompBAl_1D_spherical_2m` with its `setup_e2.py` (Phase D Al models,
products -> gas transfer below 100 bar, 28-species gas with equilibrium above
2500 K):

```
python3 setup_pontalier.py --trial 6 --out <dir> [--closure laminar|pasrd]
cd <dir> && ./Allrun && python3 summary.py
```

- The explosive phase keeps the name `compBAl`. Its products are the
  EuDetCode products-only fit of the trial
  (`utilities/thermo_fitter/eudet_pontalier/`, `pontalier_charges.py`), with
  the booster lumped into the charge; reactants Murnaghan with K0 = 5 GPa.
- Al: as in the Comp-B/Al cases, 90% of the Al is a separate phase (Beckstead
  burning with gas-phase oxidiser), 10% is in the products. Glass beads: the
  particle phase with 2500 kg/m3 and Cv 840 J/(kg K), never ignites.
- Charge radius 60.5-61.1 mm from the component densities (consistent with
  the bulb). Domain 4 m: 0.5 mm cells to 0.1 m, graded to 15 mm; run to 4 ms.
- `--closure pasrd`: PaSR with the cell-size mixing time (G3c); the default
  is laminar chemistry.

Limits: the experiment is not free air beyond the positive phase at 2 m, the
fireball jets and particle bow shocks are 3D, and the 1D case cannot mix the
products with the air (G3c). BKW underestimates the NM detonation (D 5% and
p_CJ 12% low).
