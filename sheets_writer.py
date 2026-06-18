"""
sheets_writer.py — Google Sheets integration for occultation-details
----------------------------------------------------------------------
Writes star and asteroid lookup results into the spreadsheet:
  "veritas asteroid occultations"

Four tabs are written:
  • Stellar Observed Parameters  — SIMBAD + Gaia astrometry/photometry
  • StarHorse Parameters         — VizieR I/354 percentiles
  • GSP Parameters               — Gaia Apsis best-module params + CI
  • Asteroid Parameters          — JPL SBDB orbital + physical data

Authentication:
  Requires a Google service account JSON key file. Set the path via:
    export GOOGLE_SERVICE_ACCOUNT_KEY=/path/to/key.json
  or pass --key /path/to/key.json on the CLI.

Usage:
    python sheets_writer.py star   "TYC 2924-2100-1"
    python sheets_writer.py star   "UCAC4 123-456789"
    python sheets_writer.py asteroid 433
    python sheets_writer.py asteroid Eros
    python sheets_writer.py --key /path/to/key.json star "TYC 2924-2100-1"
"""

import os
import sys
import math
import warnings
import argparse
from typing import Optional, Any

warnings.filterwarnings("ignore", category=UserWarning)

# Load .env before anything else so all os.environ reads pick it up
from dotenv import load_dotenv
load_dotenv()

import gspread
from google.oauth2.service_account import Credentials

# ── local modules ────────────────────────────────────────────────────────────
from star_lookup import (
    query_simbad, query_gaia, query_starhorse, _col, _safe_float
)
from jpl_sbdb_lookup import fetch_asteroid

# ─────────────────────────────────────────────────────────────────────────────
# Configuration — all sensitive values come from .env, never hardcoded
# ─────────────────────────────────────────────────────────────────────────────

def _require_env(name: str) -> str:
    """Return the value of an environment variable, or exit with a clear error."""
    val = os.environ.get(name, "").strip()
    if not val:
        print(f"Error: environment variable {name!r} is not set.")
        print("Copy .env.example to .env and fill in your values.")
        sys.exit(1)
    return val


SCOPES = ["https://www.googleapis.com/auth/spreadsheets"]

# Tab names exactly as they appear in Google Sheets
TAB_STELLAR   = "Stellar Observed Parameters"
TAB_STARHORSE = "StarHorse Parameters"
TAB_GSP       = "GSP Parameters"
TAB_ASTEROID  = "Asteroid Parameters"
TAB_DETECTIONS = "Occultation Detections"

# Row 1 is the header; data starts at row 2
HEADER_ROW = 1

# ─────────────────────────────────────────────────────────────────────────────
# Auth
# ─────────────────────────────────────────────────────────────────────────────

def get_client(key_path: str) -> gspread.Client:
    """Return an authenticated gspread client using a service account key."""
    creds = Credentials.from_service_account_file(key_path, scopes=SCOPES)
    return gspread.authorize(creds)


def open_sheet(client: gspread.Client) -> gspread.Spreadsheet:
    spreadsheet_id = _require_env("OCCULTATION_SPREADSHEET_ID")
    return client.open_by_key(spreadsheet_id)


# ─────────────────────────────────────────────────────────────────────────────
# Generic sheet helpers
# ─────────────────────────────────────────────────────────────────────────────

def _v(value: Any, precision: int = 6) -> Any:
    """
    Convert an astropy masked/None value to a plain Python float or empty string.
    Sheets API requires plain Python scalars — not numpy types.
    """
    if value is None:
        return ""
    try:
        import numpy as np
        if np.ma.is_masked(value):
            return ""
        f = float(value)
        return "" if math.isnan(f) or math.isinf(f) else round(f, precision)
    except (TypeError, ValueError):
        s = str(value).strip()
        return "" if s in ("--", "None", "nan", "") else s


def _pm_total(pmra, pmdec) -> Any:
    """Combined proper motion magnitude in mas/yr."""
    try:
        import numpy as np
        if np.ma.is_masked(pmra) or np.ma.is_masked(pmdec):
            return ""
        pm = math.sqrt(float(pmra)**2 + float(pmdec)**2)
        return round(pm, 4)
    except Exception:
        return ""


def find_or_append_row(ws: gspread.Worksheet, id_col: int, id_value: str) -> int:
    """
    Find the row index (1-based) where column id_col matches id_value.
    If not found, appends a new row and returns its index.
    id_col is 1-based (1 = column A).
    """
    col_values = ws.col_values(id_col)  # list of all values in that column
    id_value_str = str(id_value).strip()
    for i, cell in enumerate(col_values):
        if str(cell).strip() == id_value_str:
            return i + 1  # convert 0-based list index to 1-based row number
    # Not found — find first truly empty row after the header
    next_row = max(len(col_values) + 1, HEADER_ROW + 1)
    return next_row


def write_row(ws: gspread.Worksheet, row_num: int, values: list) -> None:
    """Write a list of values starting at column A of the given row."""
    if not values:
        return
    end_col = chr(ord('A') + len(values) - 1)
    cell_range = f"A{row_num}:{end_col}{row_num}"
    ws.update(cell_range, [values], value_input_option="USER_ENTERED")

# ─────────────────────────────────────────────────────────────────────────────
# Star data extraction
# ─────────────────────────────────────────────────────────────────────────────

def extract_star_data(star_id: str) -> dict:
    """
    Run all queries for a star and return a flat dict of every value
    needed for the three star-related sheets.
    """
    print(f"\n[LOOKUP] Querying star: {star_id!r}")

    simbad = query_simbad(star_id)
    gaia_id = simbad.get("gaia_source_id")

    gaia = {}
    starhorse = None
    if gaia_id:
        gaia = query_gaia(gaia_id)
        starhorse = query_starhorse(gaia_id)
    else:
        print("[WARN] No Gaia DR3 source_id — Gaia and StarHorse queries skipped.")

    gs  = gaia.get("gaia_source")
    ap  = gaia.get("ap")

    # ── Proper motion (combined) ──────────────────────────────────────────────
    pmra  = _col(gs, "pmra")  if gs else None
    pmdec = _col(gs, "pmdec") if gs else None
    # Propagated PM uncertainty (quadrature sum of pmra_error + pmdec_error)
    # gaia_source has pmra_error and pmdec_error but we didn't select them.
    # Use parallax_error as a proxy note — PM error not in our current query.
    # We store combined PM magnitude; uncertainty left blank (not in query).
    pm_total = _pm_total(pmra, pmdec) if gs else ""

    # ── Determine best Apsis module for GSP tab ───────────────────────────────
    # Priority: GSP-Spec (if teff present) > GSP-Phot > ESP-HS
    apsis_module = ""
    teff_val = logg_val = mh_val = radius_val = ""
    teff_lo  = logg_lo  = mh_lo  = radius_lo  = ""
    teff_hi  = logg_hi  = mh_hi  = radius_hi  = ""

    if ap is not None:
        teff_spec = _v(_col(ap, "teff_gspspec"))
        teff_phot = _v(_col(ap, "teff_gspphot"))
        teff_hs   = _v(_col(ap, "teff_esphs"))

        if teff_spec != "":
            apsis_module = "GSP-Spec"
            teff_val, teff_lo, teff_hi = (
                _v(_col(ap, "teff_gspspec")),
                _v(_col(ap, "teff_gspspec_lower")),
                _v(_col(ap, "teff_gspspec_upper")),
            )
            logg_val, logg_lo, logg_hi = (
                _v(_col(ap, "logg_gspspec")),
                _v(_col(ap, "logg_gspspec_lower")),
                _v(_col(ap, "logg_gspspec_upper")),
            )
            mh_val, mh_lo, mh_hi = (
                _v(_col(ap, "mh_gspspec")),
                _v(_col(ap, "mh_gspspec_lower")),
                _v(_col(ap, "mh_gspspec_upper")),
            )
            # FLAME radius paired with GSP-Spec → use flame radius
            radius_val, radius_lo, radius_hi = (
                _v(_col(ap, "radius_flame")),
                _v(_col(ap, "radius_flame_lower")),
                _v(_col(ap, "radius_flame_upper")),
            )
        elif teff_phot != "":
            apsis_module = f"GSP-Phot ({_v(_col(ap, 'libname_gspphot'))})"
            teff_val, teff_lo, teff_hi = (
                _v(_col(ap, "teff_gspphot")),
                _v(_col(ap, "teff_gspphot_lower")),
                _v(_col(ap, "teff_gspphot_upper")),
            )
            logg_val, logg_lo, logg_hi = (
                _v(_col(ap, "logg_gspphot")),
                _v(_col(ap, "logg_gspphot_lower")),
                _v(_col(ap, "logg_gspphot_upper")),
            )
            mh_val, mh_lo, mh_hi = (
                _v(_col(ap, "mh_gspphot")),
                _v(_col(ap, "mh_gspphot_lower")),
                _v(_col(ap, "mh_gspphot_upper")),
            )
            radius_val, radius_lo, radius_hi = (
                _v(_col(ap, "radius_gspphot")),
                _v(_col(ap, "radius_gspphot_lower")),
                _v(_col(ap, "radius_gspphot_upper")),
            )
        elif teff_hs != "":
            apsis_module = "ESP-HS"
            teff_val = teff_hs
            teff_lo  = ""
            teff_hi  = _v(_col(ap, "teff_esphs_uncertainty"))  # stored as ±σ
            logg_val = _v(_col(ap, "logg_esphs"))
            logg_hi  = _v(_col(ap, "logg_esphs_uncertainty"))
            mh_val   = ""  # ESP-HS assumes solar [M/H]
            radius_val = ""

    return {
        "star_id"         : star_id,
        "gaia_id"         : gaia_id or "",
        # Stellar Observed Parameters
        "ra"              : _v(_col(gs, "ra"))              if gs else "",
        "dec"             : _v(_col(gs, "dec"))             if gs else "",
        "parallax"        : _v(_col(gs, "parallax"))        if gs else "",
        "parallax_err"    : _v(_col(gs, "parallax_error"))  if gs else "",
        "pm_total"        : pm_total,
        "pm_err"          : "",   # not in current query — add pmra_error/pmdec_error if needed
        "b_mag"           : _v(simbad.get("flux_BT")),
        "b_mag_err"       : _v(simbad.get("flux_BT_err")),
        "v_mag"           : _v(simbad.get("flux_VT")),
        "v_mag_err"       : _v(simbad.get("flux_VT_err")),
        "g_mag"           : _v(_col(gs, "phot_g_mean_mag")) if gs else "",
        "g_mag_err"       : "",   # computed from flux_over_error if needed
        # GSP Parameters
        "apsis_module"    : apsis_module,
        "teff"            : teff_val,
        "teff_lower"      : teff_lo,
        "teff_upper"      : teff_hi,
        "logg"            : logg_val,
        "logg_lower"      : logg_lo,
        "logg_upper"      : logg_hi,
        "mh"              : mh_val,
        "mh_lower"        : mh_lo,
        "mh_upper"        : mh_hi,
        "radius"          : radius_val,
        "radius_lower"    : radius_lo,
        "radius_upper"    : radius_hi,
        # StarHorse
        "sh_teff50"       : _v(_col(starhorse, "teff50"))     if starhorse else "",
        "sh_teff16"       : _v(_col(starhorse, "teff16"))     if starhorse else "",
        "sh_teff84"       : _v(_col(starhorse, "teff84"))     if starhorse else "",
        "sh_logg50"       : _v(_col(starhorse, "logg50"))     if starhorse else "",
        "sh_logg16"       : _v(_col(starhorse, "logg16"))     if starhorse else "",
        "sh_logg84"       : _v(_col(starhorse, "logg84"))     if starhorse else "",
        "sh_mh50"         : _v(_col(starhorse, "__M_H_50"))   if starhorse else "",
        "sh_mh16"         : _v(_col(starhorse, "__M_H_16"))   if starhorse else "",
        "sh_mh84"         : _v(_col(starhorse, "__M_H_84"))   if starhorse else "",
        "sh_mass50"       : _v(_col(starhorse, "mass50"))     if starhorse else "",
        "sh_mass16"       : _v(_col(starhorse, "mass16"))     if starhorse else "",
        "sh_mass84"       : _v(_col(starhorse, "mass84"))     if starhorse else "",
        "sh_AV50"         : _v(_col(starhorse, "AV50"))       if starhorse else "",
        "sh_AV16"         : _v(_col(starhorse, "AV16"))       if starhorse else "",
        "sh_AV84"         : _v(_col(starhorse, "AV84"))       if starhorse else "",
        "sh_AG50"         : _v(_col(starhorse, "AG50"))       if starhorse else "",
        "sh_AG16"         : _v(_col(starhorse, "AG16"))       if starhorse else "",
        "sh_AG84"         : _v(_col(starhorse, "AG84"))       if starhorse else "",
        # A_B is not directly in StarHorse — placeholder for calculated value
        "sh_AB_calc"      : "",
    }

# ─────────────────────────────────────────────────────────────────────────────
# Asteroid data extraction
# ─────────────────────────────────────────────────────────────────────────────

def extract_asteroid_data(asteroid_id: str) -> dict:
    """
    Run SBDB query and return a flat dict matching the Asteroid Parameters tab.
    """
    print(f"\n[LOOKUP] Querying asteroid: {asteroid_id!r}")
    data = fetch_asteroid(asteroid_id, include_physical=True)

    if "object" not in data:
        raise ValueError(f"SBDB returned no object for {asteroid_id!r}")

    obj   = data["object"]
    orbit = data.get("orbit", {})
    phys  = {p["name"]: p for p in data.get("phys_par", [])}

    def orb(name):
        """Get value from orbital elements list by element name."""
        for el in orbit.get("elements", []):
            if el["name"] == name:
                return el.get("value"), el.get("sigma")
        return None, None

    def phys_val(name):
        p = phys.get(name, {})
        return _v(p.get("value")), _v(p.get("sigma"))

    diam_val, diam_err   = phys_val("diameter")
    period_val, per_err  = phys_val("rot_per")

    a_val, a_err   = orb("a")
    q_val, q_err   = orb("q")
    e_val, e_err   = orb("e")
    i_val, i_err   = orb("i")

    return {
        "number"     : _v(obj.get("des")),
        "name"       : str(obj.get("fullname", "")).strip(),
        "diameter"   : _v(diam_val),
        "diam_err"   : _v(diam_err),
        "period_rot" : _v(period_val),
        "period_err" : _v(per_err),
        "a"          : _v(a_val),
        "a_err"      : _v(a_err),
        "q"          : _v(q_val),
        "q_err"      : _v(q_err),
        "e"          : _v(e_val),
        "e_err"      : _v(e_err),
        "i"          : _v(i_val),
        "i_err"      : _v(i_err),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Sheet writers
# ─────────────────────────────────────────────────────────────────────────────

def write_star_to_sheets(spreadsheet: gspread.Spreadsheet, d: dict) -> None:
    """Write all star data to the three star tabs."""
    star_id = str(d["star_id"])
    gaia_id = str(d["gaia_id"])

    # ── Tab 1: Stellar Observed Parameters ───────────────────────────────────
    # Columns A–N:
    # Star ID | Gaia ID | RA | Dec | parallax | err |
    # proper motion | err | B-mag | err | V-mag | err | G-mag | err
    ws = spreadsheet.worksheet(TAB_STELLAR)
    row_num = find_or_append_row(ws, id_col=1, id_value=star_id)
    write_row(ws, row_num, [
        star_id,
        gaia_id,
        d["ra"],
        d["dec"],
        d["parallax"],
        d["parallax_err"],
        d["pm_total"],
        d["pm_err"],
        d["b_mag"],
        d["b_mag_err"],
        d["v_mag"],
        d["v_mag_err"],
        d["g_mag"],
        d["g_mag_err"],
    ])
    print(f"[SHEETS] Stellar Observed Parameters → row {row_num}")

    # ── Tab 2: StarHorse Parameters ───────────────────────────────────────────
    # Columns A–T:
    # Star ID | Gaia ID |
    # SH Temp_eff 50th | 16th | 84th |
    # SH log(g) 50th | 16th | 84th |
    # SH [m/H] 50th | 16th | 84th |
    # SH Mass 50th | 16th | 84th |
    # SH A_V 50th | 16th | 84th |
    # SH A_G 50th | 16th | 84th |
    # SH A_B (calculated)
    ws = spreadsheet.worksheet(TAB_STARHORSE)
    row_num = find_or_append_row(ws, id_col=1, id_value=star_id)
    write_row(ws, row_num, [
        star_id,
        gaia_id,
        d["sh_teff50"], d["sh_teff16"], d["sh_teff84"],
        d["sh_logg50"], d["sh_logg16"], d["sh_logg84"],
        d["sh_mh50"],   d["sh_mh16"],   d["sh_mh84"],
        d["sh_mass50"], d["sh_mass16"], d["sh_mass84"],
        d["sh_AV50"],   d["sh_AV16"],   d["sh_AV84"],
        d["sh_AG50"],   d["sh_AG16"],   d["sh_AG84"],
        d["sh_AB_calc"],
    ])
    print(f"[SHEETS] StarHorse Parameters → row {row_num}")

    # ── Tab 3: GSP Parameters ─────────────────────────────────────────────────
    # Columns A–M:
    # Star ID | Gaia ID | Apsis module |
    # Teff value | lower | upper |
    # log(g) value | lower | upper |
    # [m/H] value | lower | upper |
    # Radius value  (no CI columns in header — store value only)
    ws = spreadsheet.worksheet(TAB_GSP)
    row_num = find_or_append_row(ws, id_col=1, id_value=star_id)
    write_row(ws, row_num, [
        star_id,
        gaia_id,
        d["apsis_module"],
        d["teff"],       d["teff_lower"],   d["teff_upper"],
        d["logg"],       d["logg_lower"],   d["logg_upper"],
        d["mh"],         d["mh_lower"],     d["mh_upper"],
        d["radius"],
    ])
    print(f"[SHEETS] GSP Parameters → row {row_num}")


def write_asteroid_to_sheets(spreadsheet: gspread.Spreadsheet, d: dict) -> None:
    """Write asteroid data to the Asteroid Parameters tab."""
    # Columns A–N:
    # Asteroid Number | Asteroid Name |
    # Diameter | err | Period_rot | err |
    # a | err | q | err | e | err | i | err
    ws = spreadsheet.worksheet(TAB_ASTEROID)
    row_num = find_or_append_row(ws, id_col=1, id_value=str(d["number"]))
    write_row(ws, row_num, [
        d["number"],
        d["name"],
        d["diameter"],  d["diam_err"],
        d["period_rot"],d["period_err"],
        d["a"],         d["a_err"],
        d["q"],         d["q_err"],
        d["e"],         d["e_err"],
        d["i"],         d["i_err"],
    ])
    print(f"[SHEETS] Asteroid Parameters → row {row_num}")

# ─────────────────────────────────────────────────────────────────────────────
# Batch runner — reads Occultation Detections and processes every row
# ─────────────────────────────────────────────────────────────────────────────


def _already_populated(ws: gspread.Worksheet, id_col: int, id_value: str) -> bool:
    """
    Return True if the row for id_value already has data beyond the ID column
    (i.e., at least one non-empty cell in columns B onward on that row).
    """
    col_values = ws.col_values(id_col)
    id_value_str = str(id_value).strip()
    for i, cell in enumerate(col_values):
        if str(cell).strip() == id_value_str:
            row_num = i + 1
            # Fetch the full row to check for existing data
            row_data = ws.row_values(row_num)
            # Check whether anything beyond column A is filled
            return any(str(v).strip() for v in row_data[1:])
    return False  # not present at all → not populated


def run_batch(
    spreadsheet: gspread.Spreadsheet,
    dry_run: bool = False,
    force: bool = False,
    star_only: bool = False,
    asteroid_only: bool = False,
) -> None:
    """
    Read every row from the Occultation Detections tab (columns A and B),
    then for each unique star and asteroid:
      - query the relevant data sources
      - write results into the matching tabs
      - skip entries that already have data unless --force is passed

    Args:
        spreadsheet:   open gspread Spreadsheet object
        dry_run:       print what would be written without touching Sheets
        force:         overwrite rows that already have data
        star_only:     only process stars, skip asteroids
        asteroid_only: only process asteroids, skip stars
    """
    ws_det = spreadsheet.worksheet(TAB_DETECTIONS)
    all_rows = ws_det.get_all_values()

    if len(all_rows) < 2:
        print("[BATCH] Occultation Detections tab is empty — nothing to do.")
        return

    # Skip header row; columns are 0-indexed here
    data_rows = all_rows[1:]

    # Collect unique star IDs and asteroid names, preserving encounter order
    seen_stars     = {}   # star_id  → first row number (1-based, for reporting)
    seen_asteroids = {}   # asteroid_id → first row number

    for idx, row in enumerate(data_rows, start=2):  # row 2 = first data row
        star_id     = row[0].strip() if len(row) > 0 else ""
        asteroid_num = row[1].strip() if len(row) > 1 else ""
        if star_id and star_id not in seen_stars:
            seen_stars[star_id] = idx
        if asteroid_num and asteroid_num not in seen_asteroids:
            seen_asteroids[asteroid_num] = idx

    print(f"\n[BATCH] Found {len(seen_stars)} unique star(s) and "
          f"{len(seen_asteroids)} unique asteroid(s) in Occultation Detections.")

    # ── Process stars ─────────────────────────────────────────────────────────
    if not asteroid_only:
        ws_stellar = spreadsheet.worksheet(TAB_STELLAR)
        star_ok = star_skip = star_fail = 0

        for star_id in seen_stars:
            # Check for existing data unless force mode
            if not force and _already_populated(ws_stellar, id_col=1, id_value=star_id):
                print(f"[BATCH] SKIP star {star_id!r} — already in Stellar Observed Parameters "
                      f"(use --force to overwrite)")
                star_skip += 1
                continue

            try:
                data = extract_star_data(star_id)
                if dry_run:
                    print(f"[DRY-RUN] Would write star {star_id!r}: "
                          f"Gaia={data['gaia_id'] or 'not found'}, "
                          f"Teff={data['teff'] or '—'}, "
                          f"G={data['g_mag'] or '—'}")
                else:
                    write_star_to_sheets(spreadsheet, data)
                    print(f"[BATCH] ✓ Star {star_id!r} written.")
                star_ok += 1
            except Exception as e:
                print(f"[BATCH] ✗ Star {star_id!r} failed: {e}")
                star_fail += 1

        print(f"\n[BATCH] Stars: {star_ok} written, {star_skip} skipped, {star_fail} failed.")

    # ── Process asteroids ─────────────────────────────────────────────────────
    if not star_only:
        ws_asteroid = spreadsheet.worksheet(TAB_ASTEROID)
        ast_ok = ast_skip = ast_fail = 0

        for asteroid_num in seen_asteroids:
            # Skip check: column 1 (Asteroid Number) in Asteroid Parameters tab
            if not force and _already_populated(ws_asteroid, id_col=1, id_value=asteroid_num):
                print(f"[BATCH] SKIP asteroid {asteroid_num!r} — already in Asteroid Parameters "
                      f"(use --force to overwrite)")
                ast_skip += 1
                continue

            try:
                data = extract_asteroid_data(asteroid_num)
                if dry_run:
                    print(f"[DRY-RUN] Would write asteroid {asteroid_num!r}: "
                          f"name={data['name']}, "
                          f"diam={data['diameter'] or '—'} km")
                else:
                    write_asteroid_to_sheets(spreadsheet, data)
                    print(f"[BATCH] ✓ Asteroid {asteroid_num!r} written.")
                ast_ok += 1
            except Exception as e:
                print(f"[BATCH] ✗ Asteroid {asteroid_num!r} failed: {e}")
                ast_fail += 1

        print(f"[BATCH] Asteroids: {ast_ok} written, {ast_skip} skipped, {ast_fail} failed.")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Look up star or asteroid data and write it to Google Sheets."
    )
    parser.add_argument(
        "--key",
        default=os.environ.get("GOOGLE_SERVICE_ACCOUNT_KEY", "service_account.json"),
        help="Path to the Google service account JSON key file. "
             "Defaults to $GOOGLE_SERVICE_ACCOUNT_KEY env var, then 'service_account.json'.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Query data and print it without writing to Sheets.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_star = sub.add_parser("star", help="Look up a star by TYC/UCAC4 or any SIMBAD ID.")
    p_star.add_argument("star_id", nargs="+", help="Star identifier, e.g. TYC 2924-2100-1")

    p_ast = sub.add_parser("asteroid", help="Look up an asteroid by number, name, or designation.")
    p_ast.add_argument("asteroid_id", nargs="+", help="Asteroid identifier, e.g. 433 or Eros")

    p_batch = sub.add_parser(
        "batch",
        help="Read all star/asteroid IDs from 'Occultation Detections' and populate all tabs.",
    )
    p_batch.add_argument(
        "--force",
        action="store_true",
        help="Overwrite rows that already have data (default: skip populated rows).",
    )
    p_batch.add_argument(
        "--stars-only",
        action="store_true",
        help="Only process stars; skip asteroids.",
    )
    p_batch.add_argument(
        "--asteroids-only",
        action="store_true",
        help="Only process asteroids; skip stars.",
    )

    args = parser.parse_args()

    if args.command == "star":
        star_id = " ".join(args.star_id)
        data = extract_star_data(star_id)

        if args.dry_run:
            print("\n── Extracted data (dry run) ──")
            for k, v in data.items():
                print(f"  {k:<20s}: {v}")
            return

        if not os.path.isfile(args.key):
            print(f"Error: service account key not found at {args.key!r}")
            print("Set $GOOGLE_SERVICE_ACCOUNT_KEY or use --key /path/to/key.json")
            sys.exit(1)

        client = get_client(args.key)
        ss = open_sheet(client)
        write_star_to_sheets(ss, data)
        print(f"\n✓ Star {star_id!r} written to Google Sheets.")

    elif args.command == "asteroid":
        asteroid_id = " ".join(args.asteroid_id)
        data = extract_asteroid_data(asteroid_id)

        if args.dry_run:
            print("\n── Extracted data (dry run) ──")
            for k, v in data.items():
                print(f"  {k:<20s}: {v}")
            return

        if not os.path.isfile(args.key):
            print(f"Error: service account key not found at {args.key!r}")
            print("Set $GOOGLE_SERVICE_ACCOUNT_KEY or use --key /path/to/key.json")
            sys.exit(1)

        client = get_client(args.key)
        ss = open_sheet(client)
        write_asteroid_to_sheets(ss, data)
        print(f"\n✓ Asteroid {asteroid_id!r} written to Google Sheets.")

    elif args.command == "batch":
        if args.dry_run:
            print("[BATCH] Dry-run mode — no data will be written to Sheets.\n")
            # In dry-run we still need the spreadsheet to read the detections tab
            if not os.path.isfile(args.key):
                print(f"Error: service account key not found at {args.key!r}")
                sys.exit(1)
            client = get_client(args.key)
            ss = open_sheet(client)
            run_batch(
                ss,
                dry_run=True,
                force=args.force,
                star_only=args.stars_only,
                asteroid_only=args.asteroids_only,
            )
            return

        if not os.path.isfile(args.key):
            print(f"Error: service account key not found at {args.key!r}")
            print("Set $GOOGLE_SERVICE_ACCOUNT_KEY or use --key /path/to/key.json")
            sys.exit(1)

        client = get_client(args.key)
        ss = open_sheet(client)
        run_batch(
            ss,
            dry_run=False,
            force=args.force,
            star_only=args.stars_only,
            asteroid_only=args.asteroids_only,
        )
        print("\n[BATCH] Complete.")


if __name__ == "__main__":
    main()
