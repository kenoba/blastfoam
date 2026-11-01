# 2D axisymmetric ChargeCompBAl with a seeded contact perturbation

Plan step G3b (`report/stage5_g3_plan.md`): growth of the products/air
mixing layer from a seeded perturbation of the contact surface.

## Notes

- The case of `../ChargeCompBAl_2D_axi_2m` (G2b: full E2 physics on a 2 m
  wedge, AMR, products/air mixing layer held at 2 mm), with `endTime`
  1.3 ms, the end of the unseeded G2b run that is the reference.
- `Allrun` maps the 1D solution at 60 us as in G2b. It then seeds the
  perturbation with `perturb_contact.py` before the run. In a band of
  half-width h = 8 mm around the contact radius (0.249 m at 60 us; between
  the inward-facing secondary shock at 0.236 m and the lead shock at
  0.265 m), each cell is reset to the 1D state at the displaced radius
  r - delta(psi) cos^2(pi (r - r_c)/(2 h)). The state comes from the 1D cell
  containing that radius, so primitive and conserved fields stay
  consistent. psi is the angle from the symmetry plane; the modes
  cos(n psi) have even n, so they are symmetric about the plane and smooth
  on the axis.
- With the perturbation inside the decelerating contact, it first grows as
  Rayleigh-Taylor. The re-shock by the secondary shock (about 1 ms) adds
  Richtmyer-Meshkov growth (Balakrishnan et al. 2009, phases in
  `report/external_sources/Balakrishnan09.pdf`).
- The amplitude must be well above the grid-seeded noise. In the unseeded
  G2b run the scatter of the contact radius between rays in the window
  22.5-67.5 degrees is 1.7 mm at 100 us, 4.6 mm at 500 us and 23 mm at
  1 ms. The n = 24 and 48 amplitudes are about 0.5 mm up to 500 us.
- `contact_modes.py` evaluates the written fields. Write the cell centres
  and volumes first (`postProcess -func writeCellCentres`, `writeCellVolumes`).
  On rays in the window 22.5-67.5 degrees, away from the axis and the
  symmetry plane where the fireball is mesh-dependent, it gives:
  - the contact radius and mixing width per ray;
  - the mode amplitudes over the window;
  - the CO and H2 mass of the gas.

- With `-charge` the surface of the charge is perturbed instead, and the
  run starts from 0 with the 2D detonation.
  - `setRefinedFields` sets the charge as a sphere of radius
    0.05 m + max(delta) at 0.5 mm; the cells beyond 0.05 m + delta(psi) are
    then reset to air.
  - The fixes of `ChargeCompBAl_2D_axi_cylinder` for detonation along an
    oblique charge surface are switched on (`expansionBurnRatio 0.5`,
    products/gas heat transfer scaled by `lambda.compBAl`).
- `-t0 <time>` maps the 1D solution at another time; give the name of the 1D
  time directory, e.g. 0.0002. After 60 us the band half-width is 25 mm.
- `-fine` holds the mixing layer at 1 mm (`mixingLevel 4`) instead of 2 mm.
- `-np <n>` sets the number of processes (default 6).
- Each option is written into the case dictionaries on every run, so a run
  without it resets them.

Run with, e.g.:
```
./Allrun -t0 0.0002 --modes 24:6e-3          # option 1: seeded at 200 us
./Allrun -fine --modes 24:3e-3               # option 2: 60 us, 1 mm layer
./Allrun -charge --multimode 16 64 1e-3 --seed 1   # option 3: charge surface
```
<seed> is `--modes n:amp ...` (amplitudes in m) or
`--multimode nmin nmax rms [--seed s]` (all even n in [nmin, nmax], normal
random amplitudes scaled to the rms displacement).
