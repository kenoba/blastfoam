# EuDetCode output for the Pontalier et al. (2021) charges

Products-only fits of the gelled-nitromethane charges of Pontalier, Loiseau,
Goroshin, Zhang, Frost, "Blast enhancement from metalized explosives",
Shock Waves 31 (2021) 203-230 (Table 1), generated 2026-11-03 with the
EuDetCode build in `external/eudetcode` and BKW:

```
python3 utilities/thermo_fitter/pontalier_charges.py --run 0 5 6 7 18
```

`pontalier_charges.py` builds the composition (C-4 booster lumped in, the
glass microballoons as voids, PMMA added to a temporary copy of the EuDetCode
`Data/`), the density and the inert constraint (`Al(s):0.9` for Al as in the
Comp-B/Al fit, `Al(s):1.0` for the glass beads, none for the baseline);
`trial<N>/command.txt` is the exact command. The case entries follow from
`products_only_fit.py` on the gas thermo of the 1D case
(`tutorials/blastEulerFoam/Pontalier_NM_1D_spherical/setup_pontalier.py`
does this).

BKW underestimates nitromethane: pure NM at 1137 kg/m3 gives D = 5969 m/s
and p_CJ = 11.1 GPa (measured about 6.3 km/s and 12.5-13 GPa; THEOSTAR and
JCZ3 are no closer). BKW is kept for consistency with the Comp-B/Al and TNT
fits.

| trial | charge | rho [kg/m3] | D [m/s] | p_CJ [GPa] |
|---|---|---|---|---|
| 0 | baseline | 1127.1 | 5833.6 | 10.58 |
| 5 | 20% Al (54 um) | 1267.7 | 5642.3 | 10.56 |
| 6 | 30% Al (54 um) | 1352.2 | 5547.2 | 10.53 |
| 7 | 40% Al (54 um) | 1450.0 | 5458.2 | 10.48 |
| 18 | 30% glass (66 um) | 1337.1 | 5403.9 | 9.59 |
