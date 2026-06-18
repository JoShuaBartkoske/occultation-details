"""
Star Lookup: SIMBAD + VizieR + Gaia DR3
-----------------------------------------
Given a star ID (TYC or UCAC4), this script:

  1. Queries SIMBAD for the star's basic data and resolves its Gaia DR3 source_id,
     Tycho-2 B/V photometry, coordinates, and parallax.
  2. Queries the Gaia DR3 archive (ESA TAP) for:
       - gaiadr3.gaia_source          : G, BP, RP photometry, parallax, proper motion
       - gaiadr3.astrophysical_parameters : GSP-Phot (MARCS/PHOENIX) and GSP-Spec params,
                                            FLAME mass/radius, extinction
       - gaiadr3.astrophysical_parameters_supp : per-library GSP-Phot and GSP-Spec ANN
  3. Queries VizieR (I/354, Anders et al. 2022) for StarHorse parameters.
  4. Reports the Gaia photometry and colors needed for the Kiman et al. (2024)
     surface brightness-color relation (arxiv:2406.00229).

Usage:
    python star_lookup.py TYC 1234-567-1
    python star_lookup.py UCAC4 123-456789

Dependencies:
    pip install astroquery astropy
"""

import sys
import warnings
import textwrap
from typing import Optional

# Suppress minor astropy/astroquery deprecation noise
warnings.filterwarnings("ignore", category=UserWarning)

from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from astropy.table import Table


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def fmt(value, unit: str = "", precision: int = 4, fallback: str = "—") -> str:
    """Format a potentially-None/masked table value for display."""
    if value is None:
        return fallback
    try:
        import numpy as np
        if np.ma.is_masked(value):
            return fallback
        return f"{float(value):.{precision}f}{(' ' + unit) if unit else ''}"
    except (TypeError, ValueError):
        return str(value) if str(value) not in ("--", "None", "") else fallback


def section(title: str) -> None:
    print(f"\n{'─' * 60}")
    print(f"  {title}")
    print(f"{'─' * 60}")


def row(label: str, value: str, width: int = 34) -> None:
    print(f"  {label:<{width}}: {value}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: SIMBAD — resolve the star and grab Tycho-2 photometry + Gaia ID
# ─────────────────────────────────────────────────────────────────────────────

def query_simbad(star_id: str) -> dict:
    """
    Query SIMBAD for the star.  Returns a dict with:
        main_id, ra, dec, parallax, plx_error,
        flux_B (Tycho BT), flux_V (Tycho VT),
        gaia_source_id (string or None)

    Note: astroquery >= 0.4.8 returns flux in long format (one row per filter),
    so we pivot the table into a dict keyed by filter letter.
    """
    print(f"[SIMBAD] Querying: {star_id!r}")

    sim = Simbad()
    sim.add_votable_fields("flux", "parallax", "ids", "ra", "dec", "otype")

    result = sim.query_object(star_id)
    if result is None or len(result) == 0:
        raise ValueError(f"SIMBAD: no object found for {star_id!r}")

    # The table has one row per photometric band — grab scalar fields from row 0
    r0 = result[0]

    # Build a filter → (flux, flux_err) lookup by pivoting the long-format rows
    fluxes: dict = {}   # e.g. {'B': (11.66, 0.08), 'V': (11.37, 0.09), 'G': ...}
    for row_i in result:
        filt = str(row_i["flux.filter"]).strip()
        if filt:
            fluxes[filt] = (
                _safe_float(row_i, "flux"),
                _safe_float(row_i, "flux_err"),
            )

    def _flux(band):
        return fluxes.get(band, (None, None))[0]

    def _flux_err(band):
        return fluxes.get(band, (None, None))[1]

    # Extract Gaia DR3 source_id from the pipe-delimited IDs string
    gaia_source_id = None
    ids_str = str(r0["ids"]) if "ids" in result.colnames else ""
    for token in ids_str.split("|"):
        token = token.strip()
        if token.startswith("Gaia DR3 "):
            gaia_source_id = token.replace("Gaia DR3 ", "").strip()
            break
        if token.startswith("Gaia DR2 ") and gaia_source_id is None:
            gaia_source_id = token.replace("Gaia DR2 ", "").strip()
            print("[SIMBAD] Warning: only Gaia DR2 ID found in SIMBAD; DR3 query may fail.")

    return {
        "main_id"       : str(r0["main_id"]),
        "otype"         : str(r0["otype"]),
        "ra_deg"        : _safe_float(r0, "ra"),
        "dec_deg"       : _safe_float(r0, "dec"),
        "parallax_mas"  : _safe_float(r0, "plx_value"),
        "parallax_err"  : _safe_float(r0, "plx_err"),
        "flux_BT"       : _flux("B"),
        "flux_BT_err"   : _flux_err("B"),
        "flux_VT"       : _flux("V"),
        "flux_VT_err"   : _flux_err("V"),
        "flux_G_simbad" : _flux("G"),
        "all_fluxes"    : fluxes,         # full filter dict for reference
        "gaia_source_id": gaia_source_id,
        "all_ids"       : ids_str,
    }


def _safe_float(row, colname: str) -> Optional[float]:
    try:
        import numpy as np
        val = row[colname]
        if np.ma.is_masked(val):
            return None
        return float(val)
    except Exception:
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Gaia DR3 archive — photometry + astrophysical parameters
# ─────────────────────────────────────────────────────────────────────────────

def query_gaia(source_id: str) -> dict:
    """
    Query three Gaia DR3 tables for a given source_id (as a string integer).
    Returns a nested dict.
    """
    print(f"[GAIA]   Querying source_id = {source_id}")
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    Gaia.ROW_LIMIT = 10

    sid = int(source_id)

    # ── gaia_source: astrometry + photometry ──────────────────────────────────
    q_main = f"""
    SELECT
        source_id,
        ra, dec,
        parallax, parallax_error,
        pmra, pmdec,
        phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
        bp_rp, bp_g, g_rp,
        phot_g_mean_flux_over_error,
        phot_bp_mean_flux_over_error,
        phot_rp_mean_flux_over_error,
        ruwe,
        radial_velocity, radial_velocity_error
    FROM gaiadr3.gaia_source
    WHERE source_id = {sid}
    """

    # ── astrophysical_parameters: main Apsis table ────────────────────────────
    q_ap = f"""
    SELECT
        source_id,
        -- GSP-Phot (best-fit library: MARCS, PHOENIX, A, or OB)
        teff_gspphot, teff_gspphot_lower, teff_gspphot_upper,
        logg_gspphot, logg_gspphot_lower, logg_gspphot_upper,
        mh_gspphot, mh_gspphot_lower, mh_gspphot_upper,
        radius_gspphot, radius_gspphot_lower, radius_gspphot_upper,
        distance_gspphot, distance_gspphot_lower, distance_gspphot_upper,
        azero_gspphot, azero_gspphot_lower, azero_gspphot_upper,
        ag_gspphot, ag_gspphot_lower, ag_gspphot_upper,
        abp_gspphot, arp_gspphot,
        ebpminrp_gspphot, ebpminrp_gspphot_lower, ebpminrp_gspphot_upper,
        mg_gspphot,
        libname_gspphot,
        -- GSP-Spec (RVS spectra, FGK stars)
        teff_gspspec, teff_gspspec_lower, teff_gspspec_upper,
        logg_gspspec, logg_gspspec_lower, logg_gspspec_upper,
        mh_gspspec, mh_gspspec_lower, mh_gspspec_upper,
        alphafe_gspspec, alphafe_gspspec_lower, alphafe_gspspec_upper,
        flags_gspspec,
        -- ESP-HS (OBA hot stars only): vsini is only available here
        teff_esphs, teff_esphs_uncertainty,
        logg_esphs, logg_esphs_uncertainty,
        vsini_esphs, vsini_esphs_uncertainty,
        azero_esphs, ag_esphs,
        spectraltype_esphs, flags_esphs,
        -- FLAME (uses GSP-Phot Teff + luminosity)
        radius_flame, radius_flame_lower, radius_flame_upper,
        lum_flame, lum_flame_lower, lum_flame_upper,
        mass_flame, mass_flame_lower, mass_flame_upper,
        age_flame, age_flame_lower, age_flame_upper,
        evolstage_flame, flags_flame,
        bc_flame
    FROM gaiadr3.astrophysical_parameters
    WHERE source_id = {sid}
    """

    # ── astrophysical_parameters_supp: per-library GSP-Phot + GSP-Spec ANN ───
    q_supp = f"""
    SELECT
        source_id,
        -- MARCS-specific GSP-Phot
        teff_gspphot_marcs, logg_gspphot_marcs, mh_gspphot_marcs,
        radius_gspphot_marcs, azero_gspphot_marcs, ag_gspphot_marcs,
        -- PHOENIX-specific GSP-Phot
        teff_gspphot_phoenix, logg_gspphot_phoenix, mh_gspphot_phoenix,
        radius_gspphot_phoenix, azero_gspphot_phoenix, ag_gspphot_phoenix,
        -- OB (hot star) GSP-Phot
        teff_gspphot_ob, logg_gspphot_ob, mh_gspphot_ob,
        -- GSP-Spec ANN
        teff_gspspec_ann, logg_gspspec_ann, mh_gspspec_ann
    FROM gaiadr3.astrophysical_parameters_supp
    WHERE source_id = {sid}
    """

    results = {}
    for label, adql in [("gaia_source", q_main),
                         ("ap", q_ap),
                         ("ap_supp", q_supp)]:
        try:
            job = Gaia.launch_job(adql)
            tbl = job.get_results()
            results[label] = tbl[0] if tbl and len(tbl) > 0 else None
        except Exception as e:
            print(f"[GAIA]   Warning – {label} query failed: {e}")
            results[label] = None

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: VizieR — StarHorse (Anders et al. 2022, I/354)
# ─────────────────────────────────────────────────────────────────────────────

def query_starhorse(source_id: str) -> Optional[object]:
    """
    Query VizieR catalogue I/354 (StarHorse for Gaia EDR3) using the
    Gaia source_id.  StarHorse uses the same source_ids as Gaia EDR3/DR3.
    """
    print(f"[VizieR] Querying StarHorse (I/354) for source_id = {source_id}")

    v = Vizier(columns=["*"], row_limit=5)
    # I/354 table is "I/354/starhorse" — filter on Source column
    try:
        result = v.query_constraints(catalog="I/354/starhorse2021",
                                     Source=source_id)
        if result and len(result) > 0 and len(result[0]) > 0:
            return result[0][0]
    except Exception as e:
        print(f"[VizieR] StarHorse query error: {e}")
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Display results
# ─────────────────────────────────────────────────────────────────────────────

def display_results(star_id: str,
                    simbad: dict,
                    gaia: dict,
                    starhorse) -> None:

    gs   = gaia.get("gaia_source")
    ap   = gaia.get("ap")
    aps  = gaia.get("ap_supp")

    # ── Identity ───────────────────────────────────────────────────────────────
    print("\n" + "═" * 60)
    print(f"  Star: {star_id}")
    print(f"  SIMBAD main ID: {simbad['main_id']}")
    print("═" * 60)

    section("Identifiers")
    row("Query ID",               star_id)
    row("SIMBAD main ID",         simbad["main_id"])
    row("Object type (SIMBAD)",   simbad.get("otype", "—"))
    row("Gaia DR3 source_id",     simbad["gaia_source_id"] or "not found")

    # ── Tycho-2 photometry ────────────────────────────────────────────────────
    section("Tycho-2 Photometry (from SIMBAD)")
    row("BT magnitude",  fmt(simbad["flux_BT"],     "mag", 3))
    row("BT error",      fmt(simbad["flux_BT_err"], "mag", 3))
    row("VT magnitude",  fmt(simbad["flux_VT"],     "mag", 3))
    row("VT error",      fmt(simbad["flux_VT_err"], "mag", 3))

    # ── Gaia astrometry & photometry ──────────────────────────────────────────
    section("Gaia DR3 Astrometry & Photometry (gaia_source)")
    if gs is not None:
        row("RA (deg)",        fmt(gs["ra"],  "deg", 6))
        row("Dec (deg)",       fmt(gs["dec"], "deg", 6))
        row("Parallax",        fmt(gs["parallax"],       "mas", 4)
                               + f"  ±{fmt(gs['parallax_error'], 'mas', 4)}")
        row("PM RA (mas/yr)",  fmt(gs["pmra"],  "", 4))
        row("PM Dec (mas/yr)", fmt(gs["pmdec"], "", 4))
        row("Radial velocity", fmt(gs["radial_velocity"], "km/s", 3))
        row("RUWE",            fmt(gs["ruwe"], "", 3))
        print()
        row("G magnitude",     fmt(gs["phot_g_mean_mag"],  "mag", 4))
        row("BP magnitude",    fmt(gs["phot_bp_mean_mag"], "mag", 4))
        row("RP magnitude",    fmt(gs["phot_rp_mean_mag"], "mag", 4))
        print()
        row("BP-RP color",     fmt(gs["bp_rp"], "mag", 4))
        row("G-BP  color",     fmt(gs["bp_g"],  "mag", 4))
        row("G-RP  color",     fmt(gs["g_rp"],  "mag", 4))
    else:
        row("Status", "Gaia source query returned no data")

    # ── Kiman et al. 2024 SBCR note ───────────────────────────────────────────
    if gs is not None:
        section("Kiman et al. (2024) SBCR Inputs  [arxiv:2406.00229]")
        print("  The SBCR requires one Gaia magnitude and one Gaia color.")
        print("  Nine relations are available (3 magnitudes × 3 colors).")
        print("  Inputs collected here:\n")
        row("G  (surface brightness proxy)",  fmt(gs["phot_g_mean_mag"],  "mag", 4))
        row("BP (surface brightness proxy)",  fmt(gs["phot_bp_mean_mag"], "mag", 4))
        row("RP (surface brightness proxy)",  fmt(gs["phot_rp_mean_mag"], "mag", 4))
        row("BP-RP",                          fmt(gs["bp_rp"], "mag", 4))
        row("G-BP  (= -(BP-G))",              fmt(gs["bp_g"],  "mag", 4))
        row("G-RP",                           fmt(gs["g_rp"],  "mag", 4))

    # ── GSP-Phot ──────────────────────────────────────────────────────────────
    section("Gaia Apsis: GSP-Phot (best library)")
    if ap is not None:
        import numpy as np
        lib_raw = _col(ap, "libname_gspphot")
        lib = str(lib_raw).strip() if lib_raw is not None else ""
        if not lib or lib in ("--", "None", "") or (hasattr(lib_raw, 'mask') and lib_raw.mask):
            lib = "—"
        row("Library used", lib)
        if lib != "—":
            _print_param_with_ci(ap, "Teff",        "teff_gspphot",    "K")
            _print_param_with_ci(ap, "log g",       "logg_gspphot",    "dex")
            _print_param_with_ci(ap, "[M/H]",       "mh_gspphot",      "dex")
            _print_param_with_ci(ap, "Radius",      "radius_gspphot",  "R☉")
            _print_param_with_ci(ap, "Distance",    "distance_gspphot","pc")
            _print_param_with_ci(ap, "A_0 (541nm)", "azero_gspphot",   "mag")
            _print_param_with_ci(ap, "A_G",         "ag_gspphot",      "mag")
            _print_param_with_ci(ap, "A_BP",        "abp_gspphot",     "mag")
            _print_param_with_ci(ap, "A_RP",        "arp_gspphot",     "mag")
            _print_param_with_ci(ap, "E(BP-RP)",    "ebpminrp_gspphot","mag")
            _print_param_with_ci(ap, "M_G (abs)",   "mg_gspphot",      "mag")
        else:
            row("Note", "GSP-Phot parameters not available for this star")
    else:
        row("Status", "No GSP-Phot data returned")

    # ── GSP-Phot per library (supp) ───────────────────────────────────────────
    if aps is not None:
        section("Gaia Apsis: GSP-Phot per-library (astrophysical_parameters_supp)")
        for lib_tag, lib_label in [("marcs", "MARCS"), ("phoenix", "PHOENIX"), ("ob", "OB")]:
            teff_col = f"teff_gspphot_{lib_tag}"
            if teff_col in aps.colnames:
                teff_val = fmt(aps[teff_col], "K", 0)
                if teff_val != "—":
                    logg_val  = fmt(aps[f"logg_gspphot_{lib_tag}"],    "dex", 3)
                    mh_val    = fmt(aps[f"mh_gspphot_{lib_tag}"],      "dex", 3)
                    r_val     = fmt(aps[f"radius_gspphot_{lib_tag}"],  "R☉",  3)
                    a0_val    = fmt(aps[f"azero_gspphot_{lib_tag}"],   "mag", 3)
                    ag_val    = fmt(aps[f"ag_gspphot_{lib_tag}"],      "mag", 3)
                    print(f"\n  [{lib_label}]")
                    row("  Teff",    teff_val)
                    row("  log g",   logg_val)
                    row("  [M/H]",   mh_val)
                    row("  Radius",  r_val)
                    row("  A_0",     a0_val)
                    row("  A_G",     ag_val)

    # ── GSP-Spec ──────────────────────────────────────────────────────────────
    section("Gaia Apsis: GSP-Spec (RVS spectroscopy, FGK stars)")
    if ap is not None:
        _print_param_with_ci(ap, "Teff",      "teff_gspspec",    "K")
        _print_param_with_ci(ap, "log g",     "logg_gspspec",    "dex")
        _print_param_with_ci(ap, "[M/H]",     "mh_gspspec",      "dex")
        _print_param_with_ci(ap, "[α/Fe]",    "alphafe_gspspec", "dex")
        flags = fmt(_col(ap, "flags_gspspec"), "", 0)
        if flags != "—":
            row("Quality flags", str(ap["flags_gspspec"]))
    else:
        row("Status", "No GSP-Spec data returned")

    if aps is not None:
        teff_ann = fmt(_col(aps, "teff_gspspec_ann"), "K", 0)
        if teff_ann != "—":
            print()
            row("  Teff (ANN)",   teff_ann)
            row("  log g (ANN)",  fmt(_col(aps, "logg_gspspec_ann"), "dex", 3))
            row("  [M/H] (ANN)",  fmt(_col(aps, "mh_gspspec_ann"),   "dex", 3))

    # ── ESP-HS ────────────────────────────────────────────────────────────────
    section("Gaia Apsis: ESP-HS (hot O/B/A stars only)")
    if ap is not None:
        teff_hs = fmt(_col(ap, "teff_esphs"), "K", 0)
        if teff_hs != "—":
            row("Spectral type (ESP)",  str(ap["spectraltype_esphs"]))
            row("Teff",                 teff_hs + "  ±" + fmt(_col(ap, "teff_esphs_uncertainty"), "K", 0))
            row("log g",                fmt(_col(ap, "logg_esphs"), "dex", 3) + "  ±" + fmt(_col(ap, "logg_esphs_uncertainty"), "dex", 3))
            row("v sin i",              fmt(_col(ap, "vsini_esphs"), "km/s", 2) + "  ±" + fmt(_col(ap, "vsini_esphs_uncertainty"), "km/s", 2))
            row("A_0 (extinction)",     fmt(_col(ap, "azero_esphs"), "mag", 3) + "  ±" + fmt(_col(ap, "azero_esphs_uncertainty"), "mag", 3))
            row("A_G",                  fmt(_col(ap, "ag_esphs"), "mag", 3) + "  ±" + fmt(_col(ap, "ag_esphs_uncertainty"), "mag", 3))
            row("flags_esphs",          str(ap["flags_esphs"]))
        else:
            row("Status", "Not an ESP-HS target (star not classified O/B/A, or G > 17.65)")
    else:
        row("Status", "No Apsis data returned")

    # ── FLAME ─────────────────────────────────────────────────────────────────
    section("Gaia Apsis: FLAME (luminosity, mass, radius, age)")
    if ap is not None:
        _print_param_with_ci(ap, "Luminosity", "lum_flame",    "L☉")
        _print_param_with_ci(ap, "Radius",     "radius_flame", "R☉")
        _print_param_with_ci(ap, "Mass",       "mass_flame",   "M☉")
        _print_param_with_ci(ap, "Age",        "age_flame",    "Gyr")
        bc   = fmt(_col(ap, "bc_flame"),       "mag", 3)
        evol = fmt(_col(ap, "evolstage_flame"), "",   0)
        flgs_raw = _col(ap, "flags_flame")
        flgs = str(flgs_raw).strip() if flgs_raw is not None else "—"
        if flgs in ("--", "None", ""):
            flgs = "—"
        if evol != "—":
            row("Evol. stage (BASTI)", evol)
        if bc != "—":
            row("Bolometric correction (BC_G)", bc)
        if flgs != "—":
            row("flags_flame", flgs)
        if all(fmt(_col(ap, c), "", 0) == "—" for c in
               ["lum_flame", "radius_flame", "mass_flame", "age_flame"]):
            row("Note", "FLAME parameters not available (requires GSP-Phot Teff)")
    else:
        row("Status", "No FLAME data returned")

    # ── StarHorse ──────────────────────────────────────────────────────────────
    section("StarHorse (Anders et al. 2022, VizieR I/354)")
    if starhorse is not None:
        _sh_row(starhorse, "Teff (16th)",   "teff16",  "K")
        _sh_row(starhorse, "Teff (50th)",   "teff50",  "K")
        _sh_row(starhorse, "Teff (84th)",   "teff84",  "K")
        _sh_row(starhorse, "log g (50th)",  "logg50",  "dex")
        _sh_row(starhorse, "logg (16th)",   "logg16",  "dex")
        _sh_row(starhorse, "logg (84th)",   "logg84",  "dex")
        _sh_row(starhorse, "[M/H] (50th)",  "__M_H_50","dex")
        _sh_row(starhorse, "Distance (50th)","dist50",  "kpc")
        _sh_row(starhorse, "A_V (50th)",    "AV50",    "mag")
        _sh_row(starhorse, "Mass (50th)",   "mass50",  "M☉")
        _sh_row(starhorse, "Age (50th)",    "age50",   "Gyr")
    else:
        row("Status", "No StarHorse data found (star may not be in I/354)")

    print("\n" + "═" * 60)
    print("  Sources:")
    print("    SIMBAD  : https://simbad.u-strasbg.fr")
    print("    Gaia DR3: https://gea.esac.esa.int/archive/")
    print("    StarHorse (I/354): Anders et al. 2022, A&A 658, A91")
    print("    SBCR    : Kiman et al. 2024, arxiv:2406.00229")
    print("═" * 60)


def _print_param_with_ci(row_data, label: str, base_col: str, unit: str) -> None:
    """Print a parameter with its lower/upper confidence interval if available."""
    val    = fmt(row_data[base_col]                   if base_col in row_data.colnames else None, unit, 3)
    lo_col = base_col + "_lower"
    hi_col = base_col + "_upper"
    lo     = fmt(row_data[lo_col] if lo_col in row_data.colnames else None, "", 3)
    hi     = fmt(row_data[hi_col] if hi_col in row_data.colnames else None, "", 3)

    if val == "—":
        return  # skip missing parameters silently

    ci = f"  [{lo}, {hi}]" if (lo != "—" and hi != "—") else ""
    row(label, f"{val}{ci}")


def _col(row_data, colname: str):
    """Safely pull a column value from an astropy Row."""
    try:
        return row_data[colname]
    except (KeyError, IndexError):
        return None


def _sh_row(sh_row, label: str, colname: str, unit: str) -> None:
    """Print one StarHorse row, handling masked/missing gracefully."""
    try:
        val = fmt(sh_row[colname], unit, 3)
        if val != "—":
            row(label, val)
    except Exception:
        pass


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    # Accept "TYC 1234-567-1", "UCAC4 123-456789", or any SIMBAD-recognized ID
    star_id = " ".join(sys.argv[1:])
    print(f"\nLooking up: {star_id!r}\n")

    # 1. SIMBAD
    try:
        simbad_data = query_simbad(star_id)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    gaia_id = simbad_data["gaia_source_id"]

    # 2. Gaia
    gaia_data = {}
    if gaia_id:
        gaia_data = query_gaia(gaia_id)
    else:
        print("[GAIA]   No Gaia DR3 source_id found in SIMBAD — skipping Gaia queries.")

    # 3. StarHorse
    starhorse_row = None
    if gaia_id:
        starhorse_row = query_starhorse(gaia_id)

    # 4. Display
    display_results(star_id, simbad_data, gaia_data, starhorse_row)


if __name__ == "__main__":
    main()
