# Extra analyses

These analyses are not part of ExtractGalaxies. Each one is switched on by its
own tag (all are off by default) and writes its own file (see
[output files](outputs.md#files-of-the-extra-analyses)). Options are listed in
the [parameter reference](parameters.md#galpost-extra-analyses).

## Common conventions

* All analyses use the centre, velocity and frame of the galaxy computed by the
  standard (ExtractGalaxies) part: z is the disc axis.
* "Galaxy" quantities honour `LoadParticleList` and `SubUsePartAssignment`;
  halo-scale quantities (kinematic decomposition potential, flows, CGM,
  satellites) use all particles.
* Disc-plane profiles use annuli of cylindrical radius R from 0 to `ExtraRmax`
  (default rgal) in `ExtraBins` bins (linear, or logarithmic with `ExtraLog`),
  and only particles with |z| < `ExtraZmax` (default 5 kpc).
* Constants: G = 4.30091 × 10⁻⁶ kpc (km/s)² Msun⁻¹; 1 yr = 3.15576 × 10⁷ s.
  (The parts that reproduce ExtractGalaxies keep its own constants.)
* Solar abundances: Asplund et al. (2009), ARA&A 47, 481. For element X,
  [X/H] = log₁₀(m_X/m_H) − log₁₀(m_X/m_H)_⊙, with m_H = m − Σ (all tracked
  species, He included). For stars, element masses refer to the initial stellar
  mass (`iM` block), as in ExtractGalaxies.
* Cold and hot gas are the MUPPI phase masses (`MCLD`, `MHOT`); single-phase
  particles are cold below 10⁵ K and hot above, as in ExtractGalaxies. H₂ is
  the derived block `MMOL` (MUPPI pressure law with normalisation `Pressure`).

## Kinematic decomposition — `KinematicDecomposition`

Circularity from the binding energy (Abadi et al. 2003, ApJ 597, 21):

1. a spherically averaged potential Φ(r) is computed from all particles within
   rvir: Φ(r) = −G [M(<r)/r + Σ_{r_i > r} m_i/r_i];
2. for each radius the energy and angular momentum of a circular orbit are
   E_c(r) = v_c²/2 + Φ(r) and j_c(r) = r v_c, with v_c² = G M(<r)/r;
3. for each star E = v²/2 + Φ(r) and ε = jz / j_c(E), where j_c(E) is
   obtained by inverting E_c(r).

Stars on circular orbits in the disc plane have ε = 1, counter-rotating ones
ε = −1. Within rgal the file `kindec.*` gives the disc mass (ε >
`DiscCircularity`, default 0.7), the spheroid mass estimated as twice the
counter-rotating mass, thin (ε > 0.8) and thick/warm (0.5 < ε ≤ 0.8) discs.

Note: the "2 × counter-rotating" spheroid is identical to the naive
ExtractGalaxies decomposition (`DecompositionMethod 0`), because the sign of ε
depends only on jz; the energy-based circularity matters for the positive
thresholds. ExtractGalaxies' own circularity (`jhist.*`) uses j_c(r) at the
position of the star instead.

## Rotational support — `KappaRot`

κ_rot = Σ m (jz/R)² / Σ m v² for the stars within `KappaRotRadius` (default
30 kpc), with velocities relative to the galaxy (Sales et al. 2010, MNRAS 409,
1541); κ_co is the same restricted to co-rotating stars (jz > 0; Correa et al.
2017, MNRAS 472, L45). Values are written in the header of `kinematics.*` (and
in `scaling.*`). κ_co ≳ 0.4 is the usual threshold for disc galaxies.

## Stellar populations — `StellarPopulations`

Mass-weighted age (cosmic time now minus formation time), total metallicity
Z = Σ metals / m (He excluded), [Fe/H], [O/Fe], [Mg/Fe] in annuli; metallicity
distribution function and the [O/Fe]–[Fe/H] plane (80 × 60 bins,
[Fe/H] ∈ [−3, 1], [O/Fe] ∈ [−0.5, 1]) for the stars within rgal.

## Disc structure — `DiscStructure`

* stellar surface density Σ★(R) in annuli;
* exponential scale length R_d from a fit of ln Σ★ between `DiscFitRmin` and
  `DiscFitRmax` (defaults: R_e and 3 R_e, with R_e the projected stellar
  half-mass radius within rgal);
* projected and 3D stellar half-mass radii within rgal;
* per annulus: median |z|, rms z and an exponential scale height from a fit of
  the vertical distribution (12 bins up to 4 × median |z|).

## Kinematic profiles — `KinematicProfiles`

For stars (mass weighted) and cold gas (cold-mass weighted): mean v_φ, the
dispersions σ_R, σ_φ, σ_z (standard deviations of the cylindrical velocity
components) and V/σ = ⟨v_φ⟩ / √((σ_R² + σ_φ² + σ_z²)/3).

## Toomre stability — `ToomreQ`

* circular velocity and epicyclic frequency from the spherically enclosed mass:
  κ² = 2Ω(Ω + dv_c/dR), Ω = v_c/R;
* Q_gas = σ_gas κ / (π G Σ_gas), with Σ_gas the cold gas surface density and
  σ_gas the vertical dispersion of the cold phase (the thermal support of the
  MUPPI cold phase, T ≈ 300 K, is negligible);
* Q_★ = σ_R κ / (3.36 G Σ★);
* two-component Q from 1/Q = 1/Q_gas + 1/Q_★ (Wang & Silk 1994, ApJ 427, 759).

## Bars — `BarStrength`

A2(R) = |Σ m e^{2iφ}| / Σ m and phase φ₂ = ½ arg(Σ m e^{2iφ}) of the stars in
annuli; A2max and its radius; the bar length is the outer edge of the last
annulus, moving outwards from A2max, where A2 > A2max/2 and the phase stays
within 10° of the phase at A2max. The header also gives A2 of all stars within
the projected half-mass radius.

## Gas flows — `SphericalFlows`

Through spherical shells of radius r_k = f_k rvir (`FlowRadii`, default 0.1,
0.25, 0.5, 1) and thickness Δ = `FlowShellWidth` × r_k: rate = Σ m v_r / Δ,
separately for v_r > 0 (outflow) and v_r < 0 (inflow), for gas, cold gas, hot
gas, stars, DM and gas-phase metals. The mass loading factor is the gas outflow
rate divided by the SFR within rgal.

## Kennicutt–Schmidt relation — `KSRelation`

Per annulus: surface densities of total, cold, molecular (MMOL) and "atomic"
(cold − H₂) gas, Σ_SFR, and depletion times t = M / SFR of cold and molecular gas.

## Metallicity gradients — `MetallicityGradients`

12 + log(O/H) by number of the gas (weighted by cold mass, and by SFR) and
stellar [Fe/H] per annulus; linear fits between `GradientRmin` and
`GradientRmax` (defaults 0.5 and 2 R_e), in dex/kpc and dex/R_e.

## Dust — `DustAnalysis`

Uses the `DSTL` and `DSTS` blocks (large and small grains), summed over their
components. Surface densities, dust-to-gas ratio and dust / gas-phase metals
(`Zs` block) per annulus; totals within rgal in the header. Whether the `Zs`
metals include the dust depends on the chemistry module: check before
interpreting the dust-to-metal ratio.

## Circumgalactic medium — `HotHaloAnalysis`

Gas between `CGMRmin` × rvir (default 0.15) and rvir: total, hot and cold
masses, hot-mass-weighted temperature (temperature of single-phase hot gas,
hot-phase temperature of MUPPI particles) and metallicity, stellar mass in the
shell; baryon fraction within rvir and the cosmic value Ω_b/Ω_m (from
`OmegaBaryon` in `GadgetParamfile`).

## Black holes — `CentralBH`

All black holes (type 5) within rgal, most massive first: ID, `BHMA` mass,
`BHMD` accretion rate converted to Msun/yr, distance from the galaxy centre,
dynamical mass. The header gives M_BH/M_bulge, with the bulge mass of the
ExtractGalaxies decomposition (or of the kinematic decomposition if
`AnalyzeGalaxies` is off).

## Satellites — `Satellites`

All SubFind subhaloes within rvir of the galaxy, except the galaxy's own
subhalo, most massive in stars first, with their masses per type (low-resolution
mass separately), SubFind Vmax and half-mass radius. Needs a SubFind catalogue
in `GroupDir` (it is loaded automatically also when galaxies are selected from
FoF groups). Note: OpenGadget catalogues produced before the fix described in
[halo finding with OpenGadget](subfind_opengadget.md#known-issues) have wrong
Vmax for the largest subhaloes.

## Scaling relations — `ScalingRelations`

One line per galaxy (`scaling.*`) with ~30 global quantities (see
[outputs](outputs.md#columns-of-scaling)); with `WriteTables` the line is also
appended to `table-scaling.gal<ID>`, building time series. Quantities whose
analysis is not enabled are NaN; the basic masses, SFR, half-mass radius,
virial mass, Vmax and the low-resolution fraction are always computed.

## Maps — `Maps`

Face-on (x–y) and edge-on (x–z) surface density maps, `MapPixels` × `MapPixels`
pixels covering ±`MapSize` (default rgal), for the components in
`MapComponents` (`stars`, `cold`, `hot`, `gas`, `dm`). Each map is saved as a
PNG image (logarithmic colour scale) and as a `.npy` array in Msun pc⁻².
