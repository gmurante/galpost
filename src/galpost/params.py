# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Parameter file for galpost.

The syntax is the one of ExtractGalaxies, so old ``EGparams.*`` files are read
unchanged:

* ``Tag value`` for INT / DOUBLE / STRING parameters;
* ``Tag`` alone for logical switches (present = on);
* ``Quantities`` ... ``End`` block, one line per profile:
  ``QTY WEIGHT TYPE DENS|AVE DIFF|CUMU``;
* ``ListOfGals`` ... ``End`` block, one line per galaxy:
  ``x y z rmin ralign rmax [LOG] [COMP_CM] [POT] [TRACK_CM] [ALIGN]``;
* lines starting with ``#`` or ``%`` are comments; text after ``;`` or ``#``
  following a value is ignored.

Tag names are matched case-insensitively (the EG documentation and code
disagree on e.g. ``GadgetParamFile`` / ``GadgetParamfile``).

Every tag is described in :data:`TAGS` (type, default, section, help); the
same table generates the documentation and drives the configuration GUI.
New galpost analyses add their own tags there: they are all off by default.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional

__all__ = ["TAGS", "Quantity", "GalaxyLine", "Params", "read_params"]

LOGICAL, INT, DOUBLE, STRING = "logical", "int", "double", "string"


@dataclass(frozen=True)
class Tag:
    name: str
    kind: str
    default: object
    section: str
    help: str
    required: bool = False


def _t(name, kind, default, section, help, required=False):
    return Tag(name, kind, default, section, help, required)


_TAG_LIST = [
    # ---- files ------------------------------------------------------------
    _t("Snapshot", STRING, None, "files",
       "Snapshot file to analyse (overridden by -snap on the command line).", True),
    _t("Snapdir", STRING, "", "files", "Directory containing the snapshot."),
    _t("Outdir", STRING, "", "files", "Output directory (default: current)."),
    _t("Fileflag", STRING, "", "files", "Flag inserted in all output file names."),
    _t("GadgetParamfile", STRING, "paramfile.txt", "files",
       "GADGET parameter file of the simulation (units, periodicity)."),
    _t("MetalsFile", STRING, "metals.dat", "files", "List of tracked chemical elements."),
    _t("Galfile", STRING, "", "files", "File with an explicit list of galaxies."),
    _t("BlockTable", STRING, "", "files",
       "Optional block table describing the snapshot blocks (see gadget_io)."),
    _t("Verbose", LOGICAL, False, "files", "More information on stdout."),
    # ---- geometry of profiles ---------------------------------------------
    _t("SurfaceDensity", LOGICAL, False, "profiles", "Cylindrical bins, surface densities."),
    _t("VolumeDensity", LOGICAL, False, "profiles", "Spherical bins, volume densities."),
    _t("GridDensity", LOGICAL, False, "profiles", "Square grid on the disc plane."),
    _t("NumberOfBins", INT, None, "profiles", "Number of bins of the profiles.", True),
    _t("RotateGalaxyFrame", INT, None, "profiles",
       "0: disc in x-y plane; 1 / -1: rotate frame x->y->z / x->z->y.", True),
    _t("CmType", INT, None, "profiles",
       "Particle type defining centre and frame; outside 0-5: cold gas + stars "
       "(6: cold gas only).", True),
    _t("Pressure", DOUBLE, None, "profiles",
       "Normalisation of the pressure law for the molecular fraction.", True),
    _t("Cut", STRING, None, "profiles",
       "Vertical semi-height of the galaxy slice (number in kpc or rmin/rmax/rvir/"
       "rgal/ralign); <=0: no cut.", True),
    _t("OnlyMultiPhase", LOGICAL, False, "profiles",
       "Exclude single-phase particles from cold and hot masses."),
    _t("Htrue", LOGICAL, False, "profiles", "Convert quantities to the true h."),
    _t("DecompositionMethod", INT, None, "analysis",
       "Bulge/disc decomposition: -1 naive w/o halo, 0 naive, 1 Scannapieco, "
       "2 Murante.", True),
    # ---- analyses / outputs -------------------------------------------------
    _t("AnalyzeGalaxies", LOGICAL, False, "analysis", "Standard analysis of each galaxy."),
    _t("AlignWithAngularMomentum", LOGICAL, False, "analysis",
       "Galaxy frame from the angular momentum instead of the inertia tensor."),
    _t("AquilaAlignment", LOGICAL, False, "analysis",
       "Alignment radius = 3 x half-mass radius of CmType particles within rgal."),
    _t("WriteTables", LOGICAL, False, "output", "Append one line per snapshot to tables."),
    _t("ToFindGalaxy", LOGICAL, False, "output",
       "Write average position, velocity and rotation matrix in code units."),
    _t("WriteFullCatalog", LOGICAL, False, "output", "ASCII catalogue of all particles."),
    _t("WriteStarCatalog", LOGICAL, False, "output", "ASCII catalogue of all stars."),
    _t("WriteGalaxyCatalog", LOGICAL, False, "output",
       "ASCII catalogues of the particles of each galaxy."),
    _t("WriteGalaxySnapshot", LOGICAL, False, "output",
       "Snapshot centred on each galaxy, with all blocks of the original."),
    _t("GalaxySnapshotRadius", STRING, None, "output",
       "Radius of the galaxy snapshot (number or rmin/rmax/rgal/rvir/ralign)."),
    _t("GalaxySnapshotInCodeUnits", LOGICAL, False, "output",
       "Write the galaxy snapshot in code units."),
    _t("CubicGalaxySnapshot", LOGICAL, False, "output", "Cubic galaxy snapshot."),
    _t("AlignGalaxySnapshot", LOGICAL, False, "output",
       "Positions/velocities of the galaxy snapshot in the galaxy frame."),
    _t("LoadParticleList", STRING, "", "output", "Restrict the analysis to these IDs."),
    _t("WriteParticleList", LOGICAL, False, "output",
       "Write particle ID lists of each galaxy (also used by progenitor tracking)."),
    # ---- galaxy selection: FoF ---------------------------------------------
    _t("ReadFofHalos", LOGICAL, False, "selection", "Galaxies from FoF groups."),
    _t("FofMinMass", DOUBLE, None, "selection",
       "Minimum FoF mass (Msun); negative: number of most massive groups."),
    _t("FofRmin", DOUBLE, None, "selection", "rmin (kpc)."),
    _t("FofRmax", DOUBLE, None, "selection", "rmax (kpc); 0 or 'rvir': virial radius."),
    _t("FofRgal", DOUBLE, None, "selection", "rgal (kpc); 0: rvir/10."),
    _t("FofRalign", DOUBLE, None, "selection", "ralign (kpc); 0: rgal."),
    _t("FofRlog", LOGICAL, False, "selection", "Logarithmic radial bins."),
    _t("FofDoNotAlign", LOGICAL, False, "selection", "Do not align with the galaxy."),
    # ---- galaxy selection: SubFind -----------------------------------------
    _t("ReadSubFindGals", LOGICAL, False, "selection", "Galaxies from SubFind subhaloes."),
    _t("SubMinStellarMass", DOUBLE, None, "selection",
       "Minimum stellar mass (Msun); negative: number of most massive galaxies."),
    _t("SubRmin", DOUBLE, None, "selection", "rmin (kpc)."),
    _t("SubRmax", DOUBLE, None, "selection", "rmax (kpc); 0: R200c of the parent group."),
    _t("SubRalign", DOUBLE, None, "selection", "ralign (kpc); 0: rgal."),
    _t("SubRlog", LOGICAL, False, "selection", "Logarithmic radial bins."),
    _t("SubDoNotAlign", LOGICAL, False, "selection", "Do not align with the galaxy."),
    _t("SubUsePartAssignment", LOGICAL, False, "selection",
       "Use SubFind particle membership instead of apertures."),
    # ---- galpost additions ----------------------------------------------------
    _t("GroupDir", STRING, "", "selection",
       "Directory containing groups_NNN (default: snapshot directory)."),
    _t("TrackFile", STRING, "", "selection",
       "Track file written by galpost-track: the galaxy of each snapshot is taken from it "
       "(overrides ReadFofHalos / ReadSubFindGals / ListOfGals)."),
    _t("MaxLowResFraction", DOUBLE, -1.0, "selection",
       "Zoom runs: discard galaxies whose low-resolution mass fraction within rvir "
       "exceeds this value (<0: no check, as ExtractGalaxies)."),
    _t("LowResTypes", STRING, "2,3", "selection",
       "Particle types of low-resolution particles (comma separated)."),
    # ---- galpost extra analyses (all off by default) ------------------------
    _t("ExtraBins", INT, 30, "extra", "Number of radial bins of the extra profiles."),
    _t("ExtraRmax", STRING, "rgal", "extra",
       "Outer radius of the extra disc-plane profiles (kpc or rmin/rmax/rvir/rgal/ralign)."),
    _t("ExtraZmax", DOUBLE, 5.0, "extra",
       "Vertical semi-thickness (kpc) of the slab used by disc-plane analyses."),
    _t("ExtraLog", LOGICAL, False, "extra", "Logarithmic bins for the extra profiles."),
    _t("KinematicDecomposition", LOGICAL, False, "extra",
       "Circularity from the binding energy, eps = jz/jcirc(E) (Abadi et al. 2003), "
       "with disc/spheroid masses and histograms (files kindec.*, jhistE.*)."),
    _t("DiscCircularity", DOUBLE, 0.7, "extra", "Circularity threshold for disc stars."),
    _t("KappaRot", LOGICAL, False, "extra",
       "kappa_rot and kappa_co (Sales et al. 2010, Correa et al. 2017) of stars (in kinematics.*)."),
    _t("KappaRotRadius", DOUBLE, 30.0, "extra", "Aperture (kpc) for kappa_rot."),
    _t("StellarPopulations", LOGICAL, False, "extra",
       "Radial profiles of stellar age, metallicity, [Fe/H], [O/Fe], [Mg/Fe]; MDF and "
       "[O/Fe]-[Fe/H] histogram (files stellarpop.*, mdf.*, ofe_feh.*)."),
    _t("DiscStructure", LOGICAL, False, "extra",
       "Stellar surface density, exponential scale length, scale heights, projected "
       "half-mass radius (file disc.*)."),
    _t("DiscFitRmin", STRING, "auto", "extra",
       "Inner radius (kpc) of the exponential fit; auto: projected half-mass radius."),
    _t("DiscFitRmax", STRING, "auto", "extra",
       "Outer radius (kpc) of the exponential fit; auto: 3 x projected half-mass radius."),
    _t("KinematicProfiles", LOGICAL, False, "extra",
       "Mean rotation, sigma_R, sigma_phi, sigma_z and V/sigma of stars and cold gas "
       "(file kinematics.*)."),
    _t("ToomreQ", LOGICAL, False, "extra",
       "Toomre Q of cold gas, stars and two-component (file toomre.*)."),
    _t("BarStrength", LOGICAL, False, "extra",
       "m=2 Fourier amplitude A2(R), phase, A2max and bar length (file bar.*)."),
    _t("SphericalFlows", LOGICAL, False, "extra",
       "Inflow/outflow rates through spherical shells, per component, metal flows and "
       "mass loading (file flows.*)."),
    _t("FlowRadii", STRING, "0.1,0.25,0.5,1.0", "extra",
       "Shell radii in units of rvir (comma separated)."),
    _t("FlowShellWidth", DOUBLE, 0.1, "extra", "Shell thickness in units of its radius."),
    _t("KSRelation", LOGICAL, False, "extra",
       "Resolved Kennicutt-Schmidt relation: Sigma_SFR, Sigma_gas, Sigma_H2, depletion "
       "times (file ks.*)."),
    _t("MetallicityGradients", LOGICAL, False, "extra",
       "12+log(O/H) of cold gas and [Fe/H] of stars vs radius, with linear fits "
       "(file gradients.*)."),
    _t("GradientRmin", STRING, "auto", "extra", "Inner fit radius (kpc); auto: 0.5 R_e."),
    _t("GradientRmax", STRING, "auto", "extra", "Outer fit radius (kpc); auto: 2 R_e."),
    _t("DustAnalysis", LOGICAL, False, "extra",
       "Dust masses (DSTL/DSTS blocks), dust-to-gas and dust-to-metal profiles (file dust.*)."),
    _t("HotHaloAnalysis", LOGICAL, False, "extra",
       "Circumgalactic gas between CGMRmin*rvir and rvir: masses by phase, temperature, "
       "metallicity, baryon fraction (file cgm.*)."),
    _t("CGMRmin", DOUBLE, 0.15, "extra", "Inner radius of the CGM in units of rvir."),
    _t("CentralBH", LOGICAL, False, "extra",
       "Black holes within rgal: mass, accretion rate, offset from the centre (file bh.*)."),
    _t("ScalingRelations", LOGICAL, False, "extra",
       "One line per galaxy with global quantities for scaling relations "
       "(file scaling.*; with WriteTables also table-scaling.gal*)."),
    _t("Satellites", LOGICAL, False, "extra",
       "SubFind subhaloes within rvir of the galaxy (file satellites.*; needs a SubFind "
       "catalogue in GroupDir)."),
    _t("Maps", LOGICAL, False, "extra",
       "Face-on and edge-on surface density maps (PNG + .npy, files map_*)."),
    _t("MapSize", STRING, "rgal", "extra", "Half side of the maps (kpc or a radius name)."),
    _t("MapPixels", INT, 256, "extra", "Pixels per side of the maps."),
    _t("MapComponents", STRING, "stars,cold,hot", "extra",
       "Components of the maps: stars, cold, hot, gas, dm."),
]

TAGS: Dict[str, Tag] = {t.name: t for t in _TAG_LIST}
_LOWER = {k.lower(): k for k in TAGS}
_BLOCKS = {"quantities", "listofgals"}


@dataclass
class Quantity:
    quantity: str
    weight: str          # 'NULL' if none
    ptype: int
    dens: bool           # True: DENS, False: AVE
    cumulative: bool     # True: CUMU, False: DIFF

    @property
    def flag(self) -> int:          # EG convention: i1 + 2*i2
        return int(self.dens) + 2 * int(self.cumulative)


@dataclass
class GalaxyLine:
    x: float
    y: float
    z: float
    rmin: float
    ralign: float
    rmax: float
    log: bool = False
    compute_cm: bool = False
    pot: bool = False
    track_cm: bool = False
    align: bool = False


@dataclass
class Params:
    values: Dict[str, object] = field(default_factory=dict)
    quantities: List[Quantity] = field(default_factory=list)
    galaxies: List[GalaxyLine] = field(default_factory=list)
    path: str = ""

    def __getitem__(self, key):
        key = _LOWER.get(key.lower(), key)
        if key in self.values:
            return self.values[key]
        return TAGS[key].default

    def __contains__(self, key):
        return _LOWER.get(key.lower(), key) in self.values

    def get(self, key, default=None):
        key = _LOWER.get(key.lower(), key)
        return self.values.get(key, default if key not in TAGS else TAGS[key].default)

    def set(self, key, value):
        self.values[_LOWER.get(key.lower(), key)] = value

    def check(self):
        """Raise ValueError for missing required tags or inconsistent choices."""
        missing = [t.name for t in TAGS.values() if t.required and t.name not in self.values]
        if missing:
            raise ValueError("missing required parameters: " + ", ".join(missing))
        geo = [k for k in ("SurfaceDensity", "VolumeDensity", "GridDensity") if self[k]]
        if len(geo) != 1:
            raise ValueError("exactly one of SurfaceDensity, VolumeDensity, GridDensity "
                             f"must be given (found: {geo or 'none'})")
        if self["ReadFofHalos"]:
            for k in ("FofMinMass", "FofRmin", "FofRmax", "FofRgal", "FofRalign"):
                if k not in self.values:
                    raise ValueError(f"ReadFofHalos requires {k}")
        if self["ReadSubFindGals"]:
            for k in ("SubMinStellarMass", "SubRmin", "SubRmax", "SubRalign"):
                if k not in self.values:
                    raise ValueError(f"ReadSubFindGals requires {k}")
        if self["WriteGalaxySnapshot"] and "GalaxySnapshotRadius" not in self.values:
            raise ValueError("WriteGalaxySnapshot requires GalaxySnapshotRadius")

    def to_text(self) -> str:
        """Serialise in the paramfile syntax (used by the GUI)."""
        out = []
        for name, t in TAGS.items():
            if name not in self.values:
                continue
            v = self.values[name]
            if t.kind == LOGICAL:
                if v:
                    out.append(name)
            else:
                out.append(f"{name:<28s}{v}")
        if self.quantities:
            out.append("\nQuantities\n# Quantity Weight Type Density Cumulate")
            for q in self.quantities:
                out.append(f"{q.quantity:<4s} {q.weight:<4s} {q.ptype} "
                           f"{'DENS' if q.dens else 'AVE '} {'CUMU' if q.cumulative else 'DIFF'}")
            out.append("End")
        if self.galaxies:
            out.append("\nListOfGals")
            for g in self.galaxies:
                flags = [n for n, on in (("LOG", g.log), ("COMP_CM", g.compute_cm),
                                         ("POT", g.pot), ("TRACK_CM", g.track_cm),
                                         ("ALIGN", g.align)) if on]
                out.append(f"{g.x} {g.y} {g.z} {g.rmin} {g.ralign} {g.rmax} " + " ".join(flags))
            out.append("End")
        return "\n".join(out) + "\n"


def _convert(tag: Tag, raw: str, where: str):
    if tag.kind == INT:
        try:
            return int(raw)
        except ValueError:
            return int(float(raw))
    if tag.kind == DOUBLE:
        try:
            return float(raw)
        except ValueError:
            # 'rvir' for Fof/Sub radii means "use the default", i.e. 0 in EG
            if raw.lower() == "rvir" and tag.name in ("FofRmax", "SubRmax"):
                return 0.0
            # EG used atof(): other non-numeric values read as 0
            warnings.warn(f"{where}: {tag.name} '{raw}' is not a number, taken as 0 "
                          "(ExtractGalaxies convention)")
            return 0.0
    return raw


def _parse_quantity(line: str, where: str) -> Optional[Quantity]:
    w = line.split()
    if len(w) < 5:
        warnings.warn(f"{where}: badly formatted quantity '{line}', skipped")
        return None
    qty, wgt, typ, s1, s2 = w[:5]
    try:
        ptype = int(typ)
    except ValueError:
        warnings.warn(f"{where}: bad particle type in '{line}', skipped")
        return None
    if not 0 <= ptype <= 5:
        warnings.warn(f"{where}: wrong particle type {ptype} for {qty}, skipped")
        return None
    dens = {"DENS": True, "1": True, "AVE": False, "0": False}.get(s1)
    if dens is None:
        warnings.warn(f"{where}: wrong tag {s1}, assuming AVE")
        dens = False
    cumu = {"CUMU": True, "1": True, "DIFF": False, "0": False}.get(s2)
    if cumu is None:
        warnings.warn(f"{where}: wrong tag {s2}, assuming DIFF")
        cumu = False
    return Quantity(qty, wgt, ptype, dens, cumu)


def _parse_galaxy(line: str, where: str) -> Optional[GalaxyLine]:
    w = line.split()
    if len(w) < 6:
        warnings.warn(f"{where}: badly formatted galaxy line '{line}', skipped")
        return None
    g = GalaxyLine(*[float(x) for x in w[:6]])
    flags = {f.upper() for f in w[6:]}
    g.log, g.compute_cm, g.pot = "LOG" in flags, "COMP_CM" in flags, "POT" in flags
    g.track_cm, g.align = "TRACK_CM" in flags, "ALIGN" in flags
    return g


def _strip(line: str) -> str:
    for c in (";", "#"):
        i = line.find(c)
        if i > 0:
            line = line[:i]
    return line.strip()


def read_params(path: str, check: bool = True) -> Params:
    p = Params(path=path)
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        raw = lines[i]
        i += 1
        s = raw.strip()
        if not s or s[0] in "#%":
            continue
        s = _strip(s)
        w = s.split()
        key = w[0].lower()
        where = f"{path}:{i}"
        if key in _BLOCKS:
            while i < len(lines):
                b = lines[i].strip()
                i += 1
                if not b or b[0] in "#%":
                    continue
                if b.split()[0] == "End":
                    break
                b = _strip(b)
                if key == "quantities":
                    q = _parse_quantity(b, f"{path}:{i}")
                    if q:
                        p.quantities.append(q)
                else:
                    g = _parse_galaxy(b, f"{path}:{i}")
                    if g:
                        p.galaxies.append(g)
            continue
        if key not in _LOWER:
            warnings.warn(f"{where}: unknown tag '{w[0]}' ignored")
            continue
        tag = TAGS[_LOWER[key]]
        if tag.name in p.values:
            raise ValueError(f"{where}: tag {tag.name} given twice")
        if tag.kind == LOGICAL:
            p.values[tag.name] = True
        else:
            if len(w) < 2:
                raise ValueError(f"{where}: tag {tag.name} has no argument")
            p.values[tag.name] = _convert(tag, w[1], where)
    if check:
        p.check()
    return p


# --------------------------------------------------------------------------
_SECTION_TITLES = {
    "files": "Files and general settings",
    "profiles": "Profiles and galaxy frame",
    "analysis": "Standard galaxy analysis",
    "output": "Additional outputs",
    "selection": "Galaxy selection",
    "extra": "galpost extra analyses",
}


def markdown_reference() -> str:
    """Markdown table of all parameters (used to generate docs/parameters.md)."""
    out = ["# Parameter reference", "",
           "<!-- generated by `python -m galpost.params --markdown`; do not edit by hand -->", "",
           "Logical parameters are switched on by writing the tag alone on a line. "
           "Tags are case-insensitive. Parameters marked **required** must always be given. "
           "Besides the tags, a parameter file may contain the blocks `Quantities ... End` "
           "and `ListOfGals ... End` (see [the parameter file](parameter_file.md)).", ""]
    for sec, title in _SECTION_TITLES.items():
        tags = [t for t in _TAG_LIST if t.section == sec]
        if not tags:
            continue
        out += [f"## {title}", "", "| Tag | Type | Default | Description |",
                "|---|---|---|---|"]
        for t in tags:
            if t.kind == LOGICAL:
                default = "off"
            elif t.default is None:
                default = "**required**" if t.required else "—"
            elif t.default == "":
                default = "empty"
            else:
                default = f"`{t.default}`"
            out.append(f"| `{t.name}` | {t.kind} | {default} | {t.help} |")
        out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    import sys as _sys
    if "--markdown" in _sys.argv:
        print(markdown_reference())
