# G3: Cylindrical charges and fireball mixing (Stage 5, WP 5.2–5.4)

Plan step G3 of `status_review_2026-10-14.md`: cylindrical charges (L/D) and
mixing (RMI/RTI, LES vs RANS, PaSR). Builds on G2/G2b
(`tutorials/blastEulerFoam/ChargeCompBAl_2D_axi_2m`: full E2 physics on a
2 m wedge with AMR and the 2 mm mixing-layer criterion) and on the axis study
of 2026-10-23 (the 0 and 90 degree rays are mesh-dependent in the fireball;
grid noise near the axis is a few per cent of the radial velocity).

## G3a: Cylindrical charges (WP 5.2)

- Case: `tutorials/blastEulerFoam/ChargeCompBAl_2D_axi_cylinder`, the G2
  tutorial with a cylindrical charge of the same volume (0.953 kg) about the
  axis, set and refined to 0.5 mm in parallel after `decomposePar`
  (`setRefinedFields`), initiated at the centre (the mid-plane is then a
  symmetry plane, so the quarter plane is kept). `./Allrun <L/D>`.
- Runs: L/D = 1 and 4 first (largest contrast), then 0.5 and 2; to 1.3 ms
  (the end of the G2b sphere run, which is the reference of the same mass and
  physics).
- Diagnostics: peak overpressure and positive impulse against the angle
  (every 15 degrees) at 0.3, 0.5, 1.0 and 1.5 m; side-on (0 degrees) over
  end-on (90 degrees) ratio; ratio to the sphere at the same radius; Al cloud
  and fireball shape.
- End initiation (half plane, twice the cells) after the centre-initiated
  set, if needed for the comparison with the literature.
- Criteria: the trends of the literature for centrally initiated cylinders
  (peak and impulse side-on rise with L/D and fall end-on; side-on peaks up
  to about three times the sphere's in the near field, approaching the
  sphere at large scaled distance) and, where data for comparable scaled
  distances exist, the side-on/end-on ratios. References: Knock & Davies,
  Shock Waves 23 (2013) (spherical equivalence of cylinders), Xiao, Andrae &
  Gebbeken (effect of charge shape and initiation, free air), Anastacio &
  Knock. Values to be taken from the papers before acceptance.
- End-on values are on the axis, where the fireball is mesh-dependent (axis
  study); the air shock there agreed with 1D to a few per cent in G2b, which
  is well below the expected side/end ratios.

## G3b: Seeded mixing-layer growth (WP 5.3)

- The G2b sphere with a seeded perturbation of the products/air contact, to measure:
  - the growth of the mixing layer (Rayleigh–Taylor during the deceleration, Richtmyer–Meshkov at the re-shock);
  - its effect on the CO and Al burnout.
- Compare with mixing-layer growth reported for spherical and aluminised charges (e.g. Kuhl; Balakrishnan, Menon et al. 2009, `report/external_sources/Balakrishnan09.pdf`).
- Measured away from the axis and the symmetry plane only (22.5–67.5°).
- Case: `tutorials/blastEulerFoam/ChargeCompBAl_2D_axi_seeded` (G2b settings, to 1.3 ms). `perturb_contact.py` seeds modes cos(nψ), with n even and ψ the angle from the symmetry plane, in one of two ways:
  - `./Allrun [-t0 <time>] [-fine] <seed>`: the contact of the 1D solution mapped at t0 is displaced in a band around it. Each cell takes the state of the 1D cell at the displaced radius.
  - `./Allrun -charge <seed>`: the charge surface is perturbed at t = 0, and the run starts from 0 with the 2D detonation (with the cylinder fixes `expansionBurnRatio` and λ-scaled heat transfer).
  - `-fine` holds the mixing layer at 1 mm instead of 2 mm; `-np` sets the number of processes.
- Diagnostics (`contact_modes.py`), on rays in the window 22.5–67.5°:
  - contact radius R(ψ): the edge of the products core plus the products path length beyond it;
  - mixing width W(ψ) = ∫ 4f(1 − f) dr, with f the products fraction from the C/(C + N) element ratio of the gas;
  - mode amplitudes over the window;
  - CO and H2 mass of the gas.
- 1D contact (G2b 1D run): Rayleigh–Taylor unstable from 60 µs to 1.3 ms. Inside are denser products, the contact decelerates, and the Atwood number runs from −0.25 (60–200 µs) to −0.05 (1.3 ms). The contact stops at 0.755 m at 0.74 ms and then moves inwards.
- Unseeded G2b, in the window (the grid-seeded noise):

  | t | Scatter of R between rays | n = 24 amplitude | n = 48 amplitude |
  |---|---|---|---|
  | 100 µs | 1.7 mm | about 0.5 mm | about 0.5 mm |
  | 200 µs | 2.4 mm | 0.8 mm | 0.7 mm |
  | 500 µs | 4.6 mm | 0.5 mm | 1.1 mm |
  | 1 ms | 23 mm | 5.4 mm | 6.4 mm |
  | 1.3 ms | 34 mm | 6.0 mm | 5.8 mm |

  The rise after 0.5 ms is the mesh-seeded fingers.
- Test (2026-11-01), n = 24 with 3 mm seeded at 60 µs (band ±8 mm between the secondary shock and the lead shock):
  - The amplitude was 2.29 / 1.56 / 0.99 / 1.79 / 2.49 mm at 60 / 75 / 100 / 125 / 150 µs. The unseeded run had −0.56 mm at 100 µs and 0.07 mm at 150 µs.
  - So the amplitude first dropped to 43%, then grew again, with an e-folding time of about 55 µs. Linear Rayleigh–Taylor theory gives about 40 µs.
  - Two possible causes of the drop, not separated: a transient of the displacement-only seed (no matching velocity perturbation), or numerical diffusion while the contact moves at 2–2.6 km/s across 2 mm cells.
  - Options:
    - (1) seed later, when the contact is slower: t0 = 200 µs (contact at 0.495 m, secondary and lead shock 40 and 70 mm away), band ±25 mm, 6 mm (n = 24: 5.4 mm on the contact against 0.8 mm grid noise);
    - (2) seed at 60 µs with the mixing layer at 1 mm (`-fine`);
    - (3) seed the charge surface at t = 0 (as Balakrishnan et al.).
  - Option 1 runs here. Options 2 and 3 are set up for a larger machine.

## G3c: Turbulence and combustion closures (WP 5.4)

- RANS (k-omega SST) and LES (WALE) on the G3b case, and PaSR for the gas
  chemistry; afterburning sensitivity to the closure.

## Status

- 2026-10-23: plan; G3a case set up.
- 2026-10-23/24: the first L/D 1 and 4 runs needed four fixes: arrival times after AMR changes, heat transfer scaled by λ, expansion burn, and guards in the equilibrium and temperature solvers.
- 2026-11-01: both runs parked for continuation on a larger machine, about 10 µs per 30 min on 6 cores with 2 mm cells. The run directories and restart instructions are in `results/G3a_cylinders_parked/`; the early comparison is in `results/G3a_cylinders_early/`.
- Early comparison (2026-11-01, runs parked at 267 µs (L/D 1) and 173 µs (L/D 4), peaks relative to the 1D sphere of the same mass; the 2D sphere G2b agrees with 1D within 1.2% at 0.3 m and 4% at 0.5 m on all rays, the axis included):
  - L/D 1 at 0.3–0.5 m: end-on (axis) 2.1–2.4 times the sphere, side-on 1.4–1.6, and 0.50–0.64 at 30–60°. The shock arrives there 15–37% later than for the sphere, between the side and end waves.
  - L/D 4 at 0.3 m: end-on 2.6, side-on 1.04 (1.44 at 15° and 0.5 m), 0.52–0.59 at 45–60°.
  - With centre initiation the detonation runs along the axis into both end faces, which gives the strong end-on waves.
  - At 0.3 m L/D 4 is lower side-on than L/D 1 (see the literature check below).
- Literature check (2026-11-01), against the hydrocode studies of centre-initiated TNT cylinders by Gao et al. (2022) and Xiao, Andrae & Gebbeken (2020). Full texts are in `report/external_sources/` (`gao2022.pdf`, `xiao2020.pdf`, `knock2013.pdf`, `LangranWheeler21.pdf`).
  - Sp is the peak overpressure over that of the sphere of equal mass, at the scaled distance Z from the charge centre. The probes are at Z = 0.30 and 0.51 m/kg^1/3 at 0.3 and 0.5 m (charge mass, TNT equivalence 1).
  - L/D 4 against Gao Eq. (9)/Table 2 (valid for L/D >= 2; `utilities/cylinder_blast_reference.py`, sphere from the Henrych fit), angles from the mid-plane:

    | r [m] | 0° | 15° | 30° | 45° | 60° | 75° | 90° (axis) |
    |---|---|---|---|---|---|---|---|
    | 0.3, run | 1.04 | 1.16 | 0.94 | 0.52 | 0.59 | 0.81 | 2.55 |
    | 0.3, Gao | 0.81 | 0.84 | 0.90 | 0.64 | 0.49 | 1.08 | 2.37 |
    | 0.5, run | – | 1.44 | – | – | – | 0.65 | 2.16 |
    | 0.5, Gao | – | 1.33 | – | – | – | 1.25 | 2.45 |

  - L/D 1 (outside the fit; read from Gao Fig. 15a and Xiao Fig. 8 at Z = 0.3–0.5):

    | Direction | Run | Gao | Xiao (Z = 0.5) |
    |---|---|---|---|
    | End-on (axis) | 2.1–2.4 | 1.9–2.1 | 3.0–3.4 |
    | Side-on | 1.4–1.6 | 1.15–1.3 | ~1.45 |
    | Minimum at intermediate angles | 0.50–0.64 at 30–60° | 0.45–0.6 | – |

  - Agreement:
    - Both studies find, for centre initiation at small Z, the axial (end-on) peak the highest: for Z < Z1 = 2.68 (L/D 1) and 0.73 (L/D 4) in Gao, Eq. (5), and for Z < 2 at L/D 4 in Xiao. The detonation runs along the axis into the end faces.
    - Both find a deep minimum in the vertex region between the side and end waves.
    - The radial peak of L/D 4 is below or near the sphere at Z = 0.3. It only rises to its maximum further out, 1.5–1.6 at Z ≈ 0.75–1.2 in Gao and about 2 at Z ≈ 0.8 in Xiao.
    - The run lies within the spread between the two studies, which differ by up to 50% from each other on the axis at Z = 0.5. The largest offsets from Gao are +0.23 side-on at 0.3 m (L/D 4) and the 75° ray at 0.5 m (0.65 against 1.25). That ray is close to the axis, where the fireball is mesh-dependent (axis study).
  - The side-on L/D trend noted above (L/D 4 below L/D 1 at 0.3 m) is therefore what both studies predict at this Z, not a defect.
    - The rise of the radial peak with L/D (Knock & Davies 2013, P ∝ M (L/D)^1/3/R^3; Simoens & Lefebvre, about 1.8 for L/D 4) comes from end-initiated data, or from data at larger Z.
    - Gao's radial Sp of L/D 4 overtakes L/D 1 only from Z ≈ 0.6–0.7.
  - Caveats:
    - The studies use TNT (JWL, ideal-gas air, no afterburning).
    - Z depends on the TNT equivalence of the Comp-B/Al. With 1.3, the 0.3 m probes fall below Gao's range and the 0.5 m axis value of the fit rises to 2.76.
    - Only peaks are compared; impulse needs the complete runs.
  - Test once the runs are continued: at 1.0 m (Z ≈ 1.0), side-on L/D 4 should be 1.5–2 times the sphere and above L/D 1 (Gao about 1.2, Xiao about 1.3). The axis peak of L/D 4 should drop below the sphere by Z ≈ 1.5–2.
