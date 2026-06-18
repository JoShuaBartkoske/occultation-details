"""
JPL Small Body Database (SBDB) Lookup
--------------------------------------
Fetches asteroid data from the JPL SBDB API for a given asteroid ID.

Usage:
    python jpl_sbdb_lookup.py <asteroid_id>

Asteroid ID can be:
    - IAU number           e.g. 433 (Eros), 1 (Ceres)
    - Provisional designation  e.g. "2015 AB"
    - Name                     e.g. Eros
    - SPK-ID                   e.g. 2000433
"""

import sys
import json
import urllib.request
import urllib.parse
import urllib.error

SBDB_API_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"


def fetch_asteroid(asteroid_id: str, include_physical: bool = True) -> dict:
    """
    Query the JPL SBDB API for a single asteroid.

    Parameters
    ----------
    asteroid_id : str
        Asteroid designation, name, IAU number, or SPK-ID.
    include_physical : bool
        If True, request physical parameters (diameter, H, albedo, etc.).

    Returns
    -------
    dict
        Parsed JSON response from the API.
    """
    params = {
        "sstr": asteroid_id,
        "phys-par": "1" if include_physical else "0",
        "alt-des": "1",
    }
    query = urllib.parse.urlencode(params)
    url = f"{SBDB_API_URL}?{query}"

    try:
        with urllib.request.urlopen(url, timeout=15) as response:
            data = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        print(f"HTTP error {e.code}: {e.reason}")
        sys.exit(1)
    except urllib.error.URLError as e:
        print(f"Network error: {e.reason}")
        sys.exit(1)

    return data


def display_object_info(data: dict, asteroid_id: str = "") -> None:
    """Print a formatted summary of the asteroid data."""

    # Handle multiple-match response
    if data.get("code") == 300:
        print("Multiple objects matched. Please refine your search term.")
        print("Matching objects:")
        for item in data.get("list", []):
            print(f"  {item.get('pdes', '?'):20s}  {item.get('name', '')}")
        return

    # Handle not-found or error
    if "object" not in data:
        print("No object data returned. Check the ID and try again.")
        print(json.dumps(data, indent=2))
        return

    obj = data["object"]
    orbit = data.get("orbit", {})
    phys = data.get("phys_par", [])

    # ── Object identity ────────────────────────────────────────────────────────
    print("=" * 60)
    print(f"  {obj.get('fullname', 'Unknown')}")
    print("=" * 60)

    kind_map = {"an": "Numbered asteroid", "au": "Unnumbered asteroid",
                "cn": "Numbered comet",   "cu": "Unnumbered comet"}
    kind = kind_map.get(obj.get("kind", ""), obj.get("kind", ""))

    print(f"  Type          : {kind}")
    print(f"  Designation   : {obj.get('des', '—')}")
    print(f"  SPK-ID        : {obj.get('spkid', '—')}")

    orbit_class = obj.get("orbit_class", {})
    print(f"  Orbit class   : {orbit_class.get('name', '—')} ({orbit_class.get('code', '—')})")
    print(f"  NEO           : {'Yes' if obj.get('neo') else 'No'}")
    print(f"  PHA           : {'Yes' if obj.get('pha') else 'No'}")

    # Alternate designations
    alt_des = obj.get("des_alt", [])
    if alt_des:
        alt_strs = [list(d.values())[0] for d in alt_des[:5]]
        print(f"  Alt. desig.   : {', '.join(alt_strs)}")

    # ── Orbital elements ───────────────────────────────────────────────────────
    if orbit:
        elements = {el["name"]: el for el in orbit.get("elements", [])}
        print()
        print("  Orbital Elements")
        print("  ----------------")

        def fmt_el(name, label=None):
            el = elements.get(name)
            if el:
                units = f" {el['units']}" if el.get("units") else ""
                sigma = f"  ±{el['sigma']}" if el.get("sigma") else ""
                lbl = label or el.get("title", name)
                print(f"  {lbl:<30s}: {el['value']}{units}{sigma}")

        fmt_el("a",  "Semi-major axis (a)")
        fmt_el("e",  "Eccentricity (e)")
        fmt_el("i",  "Inclination (i)")
        fmt_el("om", "Long. ascending node (Ω)")
        fmt_el("w",  "Arg. of perihelion (ω)")
        fmt_el("q",  "Perihelion distance (q)")
        fmt_el("ad", "Aphelion distance (Q)")
        fmt_el("per", "Orbital period")

        print()
        print(f"  {'Orbit condition code':<30s}: {orbit.get('condition_code', '—')}")
        print(f"  {'MOID (Earth)':<30s}: {orbit.get('moid', '—')} au")
        print(f"  {'First observation':<30s}: {orbit.get('first_obs', '—')}")
        print(f"  {'Last observation':<30s}: {orbit.get('last_obs', '—')}")
        print(f"  {'# observations used':<30s}: {orbit.get('n_obs_used', '—')}")
        print(f"  {'Data arc':<30s}: {orbit.get('data_arc', '—')} days")
        print(f"  {'Orbit solution date':<30s}: {orbit.get('soln_date', '—')}")

    # ── Physical parameters ────────────────────────────────────────────────────
    if phys:
        print()
        print("  Physical Parameters")
        print("  -------------------")
        for p in phys:
            name  = p.get("title") or p.get("name", "?")
            value = p.get("value", "—")
            sigma = f"  ±{p['sigma']}" if p.get("sigma") else ""
            units = f" {p['units']}" if p.get("units") else ""
            ref   = f"  [{p['ref']}]" if p.get("ref") else ""
            print(f"  {name:<30s}: {value}{units}{sigma}{ref}")

    print("=" * 60)
    print(f"  Source: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr={urllib.parse.quote(str(obj.get('des', asteroid_id)))}")
    print("=" * 60)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    asteroid_id = " ".join(sys.argv[1:])  # allow multi-word designations
    print(f"Looking up: {asteroid_id!r}")
    print()

    data = fetch_asteroid(asteroid_id)
    display_object_info(data, asteroid_id)


if __name__ == "__main__":
    main()
