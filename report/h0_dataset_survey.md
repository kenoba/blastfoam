# H0: Experimental datasets for the aluminised-charge validation (survey 2026-11-03, papers read)

Plan step H0 (`status_review_2026-10-14.md`): identify an experimental free-air blast dataset for an aluminised charge. It is needed to decide between the laminar and the PaSR limit of the afterburning (G3c: about 5% in impulse and 1–2% in peak at 0.3–0.5 m in 1D).

## What the dataset needs

- Charge composition with the Al fraction and particle size; mass, shape, density and initiation.
- Free air, or a height of burst at which the ground reflection arrives after the positive phase at the gauges. Otherwise the geometry must be modelled.
- Side-on pressure histories (peak, arrival, positive impulse) at several distances. The range Z ≈ 0.5–5 m/kg^1/3 is where the afterburning acts.
- Ideally a reference shot without Al (or with inert particles) from the same series, so that the Al contribution is measured directly.
- Ideally fireball size and temperature (the afterburning limits differ in the fireball first).

## Candidates

Papers 1, 2, 5, 6, 7 and 8 are in `external_sources/data_sources/` (2026-11-03); the DTIC reports (3) are not yet available.

| # | Source | Explosive | Arrangement | Measured | Fit for H | Access |
|---|---|---|---|---|---|---|
| 1 | Pontalier, Loiseau, Goroshin, Zhang, Frost, "Blast enhancement from metalized explosives", Shock Waves 31 (2021) 203–230, doi:10.1007/s00193-021-00994-z | Gelled NM (95.6 NM / 4 PMMA / 0.4 GMB) in a 12.7 cm glass sphere (about 1 L, 0.98–1.41 kg), with Al (13, 54, 114 µm; 10–50%), Mg (240 µm) or glass beads (66 µm; 10, 30, 50%); central C-4 booster of 86–93 g | HOB 1.52 m, free field until the ground reflection (after the positive phase at 2 m); 20 single shots, altitude 700 m (p0 ≈ 93 kPa) | Table 3: arrival times at 1.0, 1.5, 2.0 m (gauges); Table 5: gauge peak at 2.0 m; Table 7: video peak at 2.0 m; Figs. 11, 19, 20: fireball radius and video peak to about 5 m. Impulse not analysed. | **Best for the trend with Al loading** and the inert control. Scatter is large: one shot per configuration, gauge and video peaks at 2 m differ by up to 42%, and 40% Al gives 1.06–1.58 times the baseline peak depending on the particle size. Needs a gelled-NM fit and a C-4 booster. | Local copy |
| 2 | Trzciński, Cudziło, Paszula, PEP 32 (2007) 502–508, doi:10.1002/prep.200700202 | RDXph (RDX/wax 94/6) with 15, 30, 45, 60% Al (powder or flake), 30% LiF (inert), and TNT; 268 g cylinders, 50 mm diameter, 20 g HMX booster | Inside a concrete bunker, charge 1.5 m above the ground; one side-on gauge at 2 m and one reflected gauge at 2.85 m; three shots each | Peak and impulse at 2 m in figures only (Figs. 3, 4, with ranges); light duration; QSP in 0.15 and 7 m³ chambers | Closest base explosive, with RDXph, LiF and TNT references. Free-field data are thin: one distance, values to be read from figures. Result: peak 5–17% below RDXph, impulse 8–15% above. The chamber QSP is a cleaner integral test of the total afterburning, but needs a confined case. | Local copy |
| 3 | TNT equivalences of H-6, HBX-1, HBX-3, Tritonal (Swisdak, NSWC/WOL TR 75-116, DTIC ADA018544; UFC 3-340-02 Table 2-1) | H-6 (RDX/TNT/Al/wax about 45/30/20/5) | Free air, compiled | Equivalent weight for peak and impulse | Integrated check, the same way as E1. Values not yet verified. | DTIC (site down 2026-11-03) |
| 4 | Bian et al., Energ. Mater. Front. (2026), `external_sources/Wen-xiang26.pdf` | CL-20/Al/KP, 36% Al, 1 kg | 1 m HOB, 2–8 m | Peak, impulse; fireball | Low: oxidiser, ground, inconsistent data | Local copy |
| 5 | Chen, Xu, Wu, Liu, CEJEM 13 (2016) 117–134 (`data_sources/Chen.pdf`) | RDX/Al/AP, RDX/Al/B/AP, 500 g cylinders | 1.5 m HOB, 3–7 m | Peak, TNT equivalence, fireball temperature | Low: AP, small cylinders, ground | Local copy |
| 6 | Park et al., J. Korea Inst. Mil. Sci. Technol. 8(4) (2005) | Four high-Al castables, no compositions | 4–50 m | TNT equivalences only | Low | Local copy (scanned) |
| 7 | Farrimond et al., Int. J. Prot. Struct. (2024) | PE4, PE8, PE10 hemispheres (16, 7, 8 repeats) | Surface burst, far field | Pressure histories; shot-to-shot variability | No Al. **Sets the measurement noise floor:** well-controlled repeats scatter by ±2.5% in arrival and ±6–8% in peak and impulse. | Local copy |
| 8 | Zhang et al., Adv. Mater. Sci. Eng. (2022) 1751182 | RDX/Al/EVA with 20–40% Al (4–5 µm), 5% binder | 1.5 m HOB, "wall" sensors at 1.2–4.5 m | Table 2: peak, impulse, duration | Low: bare-charge mass and shape not given, no Al-free shot, impulse units garbled (labelled MPa·s), sensor mounting unclear | Local copy |

## Assessment

- **The laminar/PaSR choice cannot be decided from free-field blast data.** In 1D the two limits differ by 1–2% in peak and about 5% in impulse at 0.3–0.5 m (G3c). Repeat shots under good control scatter by ±6–8% in peak and impulse (7). The aluminised series are noisier still: single shots and up to 42% between two methods at the same distance (1). The 2D run (G3c, PaSR with the cell-size mixing time) will show whether the difference grows in 2D; unless it reaches well over 10%, the blast data only show that both limits are inside the scatter.
- **What the data can test:** the Al effect relative to an Al-free or inert-particle reference, which is 10–60% in peak and 10–30% in arrival time:
  - 2: peak 5–17% below and impulse 8–15% above RDXph at 2 m (268 g, Z ≈ 3.1 m/kg^1/3);
  - 1: inert glass delays the arrival at 2 m by up to 23% and lowers the peak; Al advances the arrival by 10–20%.
  The model must reproduce the sign and the size of these changes with active Al and with inert Al (the reacting diameter model switched off).
- **Arrival times** are the most robust measured quantity (±2.5%), and are directly sensitive to early energy release. Pontalier's Table 3 has them at 1.0, 1.5 and 2.0 m for every shot.
- The fireball, where the two limits differ most (CO burnt 7 g against about 0 in 1D), is only measured as radius (1, Fig. 11) and light duration (2); neither is a temperature or composition.

## Recommended path

1. **H1 = Pontalier (1).** Spherical, central initiation, free field, about 1 kg: it fits the existing 1D spherical setup. Needs a EuDetCode fit of the gelled NM (CHNO; the PMMA and GMB as inert or fitted), the C-4 booster lumped into the charge or modelled as a central sphere, and the Al loadings 0, 10, 30 and 40% (54 µm) plus 30% glass. Compare arrival times at 1.0–2.0 m and peaks at 2 m against Tables 3, 5 and 7, with both the laminar and the PaSR chemistry.
2. **Secondary check = Trzciński (2)** in free field: RDXph and RDXph/Al 70/30 (closest to our charge), peak and impulse at 2 m read from Figs. 3 and 4 with their ranges. The bunker, the 1.5 m height and the cylinder shape make this a 2D case or a rough 1D estimate.
3. **H-6 equivalences (3)** when DTIC is back: an integral check of the Comp-B/Al model itself, without a new fit.
4. **FOI in-house data**, if any exist, would still be better than all of these.
