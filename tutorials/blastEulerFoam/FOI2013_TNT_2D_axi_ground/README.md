# FOI June 2013 shots: 2D axisymmetric case with the ground

Second stage of the FOI 2013 comparison (`../FOI2013_TNT_1D_spherical`,
`report/h0_dataset_survey.md` candidate 9): the 1D free-air solution is mapped
at `t0` (default 1.5 ms, before the shock reaches the ground at HOB 2 m) onto a
2D axisymmetric mesh with the ground, and run until the positive phase has
passed the gauges P1-P3 (2, 4, 6 m along the ground; P3 arrival about 12 ms).

## Running

Whole chain for one shot (1D, then 2D), or a list of shots:

```
./run_foi2013.sh -shot 4 [-closure laminar|pasrd] [-np 6] [-mixingLevel 3] [-dir foi2013_runs]
./run_campaign.sh [-jobs 2] [-np 6] [shot:closure ...]
```

Both need OpenFOAM-12 and this blastfoam build sourced, and Python 3 with
numpy; the EuDetCode fits of all shots are in
`utilities/thermo_fitter/eudet_foi2013/`, so EuDetCode itself is not needed.
Run them detached (`setsid nohup ... &`): a session that ends kills its
child processes.

- `run_foi2013.sh` works in `<dir>/shot<N>_<closure>/{1D,2D}`. The 1D stage
  (serial) runs to the mapping time t0, so that t0 is written exactly, then
  continues to 5 ms (`-end1D`): the 2D refinement table uses the 1D arrival
  times at all probes. The 2D stage maps at t0 and runs to 18 ms on `-np`
  processes. A finished stage is skipped on a rerun; a stage that stopped
  part way starts again from scratch.
- `run_campaign.sh` runs `run_foi2013.sh` for several shots, `-jobs` at a
  time. Default order: HOB 2 m (shots 4 TNT, 5 TNT/Al), 0.5 m (2 and 9 TNT
  repeats, 7 TNT/Al), 1.0 m (3, 6), 0.15 m (1, 8); laminar first, then PaSR.
  For a 2D stage alone: `setup_foi2013_2d.py --oneD <1D case> --out <dir>`,
  then `./Allrun` in it.
- Mapping time (`T0_BY_HOB` in `setup_foi2013_2d.py`): 15 us, 0.1, 0.4 and
  1.5 ms for HOB 0.15, 0.5, 1.0 and 2 m (1D shock at 0.75-0.9 of the HOB);
  `setup_foi2013_2d.py` stops if the 1D shock has reached the HOB by then.
- Cost (estimate from the 0.9 ms mapping test, about 75 000 cells): roughly
  a day on 6 cores for HOB 2 m with the mixing layer at 4 mm (`-mixingLevel
  2`, 8 mm: about half). The low HOBs map earlier on finer cells (0.5 mm
  at 15 us for HOB 0.15 m) and cost more.

## Case

- Mesh: 5-degree wedge about the vertical axis, 7.04 x 5.12 m, 32 mm base
  cells; the ground is a symmetry plane (rigid, inviscid), the side and top
  are transmissive. Charge centre on the axis at the HOB (`charge.txt` of
  the 1D case).
- Physics: `constant/` of the 1D case. Numerics, initial-field set and AMR of
  `../ChargeCompBAl_2D_axi_2m` (G2b): shock and contact (20% jumps of rho,
  and absolute jumps of p above `--pressureDelta`, 2 kPa per face, coarsened
  below 1 kPa; the 20% jumps of p alone lose the weak shock beyond about
  2.5 m, which then runs on the 32 mm base cells) refined to about the 1D
  cell size at the shock radius (a time
  table from the 1D arrival times), never below `--minLevel` (2: 8 mm);
  the products/air mixing layer at `--mixingLevel` (3:
  4 mm), marked by the gas carbon/(carbon + nitrogen) element ratio (0.667
  in the TNT products, computed from the 1D transfer composition).
- `Allrun` maps the 1D solution in parallel after `decomposePar` with
  `rotateFields -refine -extend` centred at the HOB, then runs on
  `--np` processes with load balancing (`beginTime 0`: restart).
- Output: pressure at the ground every 0.5 m (P1-P3 at 2, 4, 6 m) and on rays
  up, sideways and down from the charge (1D radii); sample lines along the
  ground, the axis and the horizontal through the charge; peak-pressure field.

Test (2026-11-03, shot 4 laminar mapped at 0.9 ms on 2 processes): 35 200
base cells refined to about 75 000 at the mapping; the solver starts cleanly.

Shot 4 laminar (2026-11-13, mapped at 1.5 ms, 18 ms on 6 processes, about
45 h; 100 000 cells early, 48 000 at the end): with `--pressureDelta 2000`
the free-field peaks at 2.5 and 3 m are within 3% of the 1D run on 8 mm
cells (20% jumps only: 6-15% low). Ground gauges P1/P2/P3 (2/4/6 m): peak
200.5/74.8/39.5 kPa (measured 440.9 (doubtful)/99.9/45.1; 20% jumps only
164.4/64.6/34.8), positive impulse 122.4/72.2/48.5 Pa s (96.9/77.8/57.3),
arrival 3.664/7.528/12.254 ms (3.364/7.287/12.054).
