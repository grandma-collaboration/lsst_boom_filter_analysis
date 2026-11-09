import requests
from datetime import datetime, timedelta, timezone
from astropy.time import Time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import math
import os
import json
import csv
import warnings
import pandas as pd
# host_association is imported in main() only when host association runs (not with SKIP_HOST_ASSOCIATION),
# so SKIP_HOST_ASSOCIATION also works in envs where it cannot be imported (e.g. pandas 3).
if hasattr(pd.errors, "SettingWithCopyWarning"):  # removed in pandas 3
    warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
warnings.filterwarnings("ignore", "Mean of empty slice", category=RuntimeWarning)
try:
    pd.options.mode.chained_assignment = None
except Exception:  # option removed in pandas 3
    pass


#Cosmology when NED has no DistMpc --> to fall back to redshift
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
LOWZ_UNCERTAIN = 0.01
MAX_HOST_ARCSEC = 10.0

#galactic extinction coefficients R=A_band/E(B-V) - Schlafly & Finkbeiner (2011) ; R_V = 3.1 
EXTINCTION_R = {
    "u": 4.239, "g": 3.303, "r": 2.285, "i": 1.698, "z": 1.263, "y": 1.08,
}
FID_TO_BAND = {0: "u", 1: "g", 2: "r", 3: "i", 4: "z", 5: "y"}

def get_band(cand):
    b = cand.get("band")
    if isinstance(b, str) and b:
        return b.lower()
    fid = cand.get("fid")
    if fid is not None:
        return FID_TO_BAND.get(fid)
    return None

def galactic_extinction(ebv, band):
    if ebv is None or band is None:
        return None
    R = EXTINCTION_R.get(band)
    if R is None:
        return None
    return R * ebv

CSV_COLUMNS = [
    ("_id", "candid"), ("objectId", "objectid"), ("jd", "jd"),
    ("band", "band"), ("magpsf", "magpsf"), ("snr_psf", "snr_psf"),
    ("has_nearby_host", "has_nearby_host"),
    ("host_source", "host_source"),
    ("host_z", "host_z"), ("z_unc", "z_unc"), ("z_tech", "z_tech"),
    ("z_qual", "z_qual"),
    ("host_sep_arcsec", "host_sep_arcsec"), ("host_sep_kpc", "host_sep_kpc"),
    ("ebv", "ebv"), ("a_band", "a_band"),
    ("dm", "dm"), ("dm_source", "dm_source"),
    ("abs_mag", "abs_mag"), ("abs_mag_corr", "abs_mag_corr"),
    ("lowz_uncertain", "lowz_uncertain"),
    #GHOST columns
    ("ghost_status", "ghost_status"),
    ("ghost_ra", "ghost_ra"), ("ghost_dec", "ghost_dec"),
    ("ghost_z", "ghost_z"), ("ghost_sep_arcsec", "ghost_sep_arcsec"),
    ("ghost_objid", "ghost_objid"), ("ghost_type", "ghost_type"),
    ("ghost_error", "ghost_error"),
    #NED-LVS DLR columns
    ("nedlvs_status", "nedlvs_status"),
    ("nedlvs_ra", "nedlvs_ra"), ("nedlvs_dec", "nedlvs_dec"),
    ("nedlvs_dlr", "nedlvs_dlr"), ("nedlvs_sep_arcsec", "nedlvs_sep_arcsec"),
    ("nedlvs_objname", "nedlvs_objname"), ("nedlvs_objtype", "nedlvs_objtype"),
    ("nedlvs_z", "nedlvs_z"), ("nedlvs_distmpc", "nedlvs_distmpc"),
    ("nedlvs_n_candidates", "nedlvs_n_candidates"),
    ("nedlvs_n_ranked", "nedlvs_n_ranked"),
    ("nedlvs_diam_arcsec", "nedlvs_diam_arcsec"),
    ("nedlvs_a_arcsec", "nedlvs_a_arcsec"),
    ("nedlvs_b_arcsec", "nedlvs_b_arcsec"),
    ("nedlvs_axis_ratio", "nedlvs_axis_ratio"),
    ("nedlvs_pa_deg", "nedlvs_pa_deg"),
    #Agreement (all 3 pairwise separations)
    ("host_agreement", "host_agreement"),
    ("ghost_nedlvs_sep_arcsec",  "ghost_nedlvs_sep_arcsec"),
    ("ghost_legacy_sep_arcsec",  "ghost_legacy_sep_arcsec"),
    ("nedlvs_legacy_sep_arcsec", "nedlvs_legacy_sep_arcsec"),
    #Legacy DR10 TAP DLR
    ("legacy_status", "legacy_status"),
    ("legacy_ra", "legacy_ra"), ("legacy_dec", "legacy_dec"),
    ("legacy_dlr", "legacy_dlr"), ("legacy_sep_arcsec", "legacy_sep_arcsec"),
    ("legacy_type", "legacy_type"),
    ("legacy_a_arcsec", "legacy_a_arcsec"), ("legacy_b_arcsec", "legacy_b_arcsec"),
    ("legacy_n_candidates", "legacy_n_candidates"), ("legacy_n_ranked", "legacy_n_ranked"),
]


def _clean_for_output(value):
    try:
        import numpy as _np
        import pandas as _pd
        if value is None:
            return None
        try:
            if _pd.isna(value):
                return None
        except Exception:
            pass
        if isinstance(value, (_np.integer,)):
            return int(value)
        if isinstance(value, (_np.floating,)):
            value = float(value)
            return value if math.isfinite(value) else None
        if isinstance(value, (_np.bool_,)):
            return bool(value)
    except Exception:
        pass
    return value

def _json_default(value):
    cleaned = _clean_for_output(value)
    if cleaned is not value:
        return cleaned
    return str(value)

def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([h for _, h in CSV_COLUMNS])
        for r in rows:
            writer.writerow([_clean_for_output(r.get(k)) for k, _ in CSV_COLUMNS])

"""
Every parameter is set in scripts/ghost/ghost_inputs.py. Edit it, then run
from the root of the repository:

    python3 scripts/ghost/boom_filter_ghost.py
"""
# CONFIG
BOOM_URL = "https://api.kaboom.caltech.edu"

# AUTH
def get_token():
    username = os.environ.get("BOOM_USERNAME")
    if not username:
        raise RuntimeError("BOOM_USERNAME environment variable not set")

    password = os.environ.get("BOOM_PASSWORD")
    if not password:
        raise RuntimeError("BOOM_PASSWORD environment variable not set")

    url = f"{BOOM_URL}/auth"
    response = requests.post(
        url,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data={
            "username": username,
            "password": password,
        },
    )
    response.raise_for_status()
    data = response.json()
    token = data["access_token"]
    expires_at = datetime.now(timezone.utc) + timedelta(seconds=data.get("expires_in", 3600))
    return token, expires_at


# LSST PIPELINE
# The filter is read from FILTER_FILE in ghost_inputs.py (a JSON list, same
# format as configs/*.json), as in scripts/boom_filter_per_run.py.
def load_pipeline(filter_file):
    with open(filter_file, "r") as f:
        pipeline = json.load(f)

    if not isinstance(pipeline, list):
        raise ValueError("Filter file must contain a JSON list (Mongo pipeline)")

    return pipeline


# BOOM QUERIES
def run_filter(token, pipeline, start_time_iso, end_time_iso):
    url = f"{BOOM_URL}/filters/test"

    start_jd = Time(start_time_iso, scale='utc').jd
    end_jd = Time(end_time_iso, scale='utc').jd

    payload = {
        "pipeline": pipeline,
        "start_jd": start_jd,
        "end_jd": end_jd,
        "survey": "LSST",
        "permissions": {"LSST": []},
        "limit": 10000,
    }

    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }

    response = requests.post(url, json=payload, headers=headers)

    if response.status_code != 200:
        print(response.status_code)
        print(response.text)
        response.raise_for_status()

    return response.json()

def get_candidates(token, candid_list, batch_size=1000):
    url = f"{BOOM_URL}/queries/find"
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }
    all_results = []
    #loop over batches
    for i in range(0, len(candid_list), batch_size):
        batch_ids = candid_list[i:i + batch_size]

        payload = {
            "catalog_name": "LSST_alerts",
            "filter": {
                "_id": {"$in": batch_ids}
            },
        }

        print(f"Querying candidate batch {i // batch_size + 1} ({len(batch_ids)} ids)")

        response = requests.post(url, json=payload, headers=headers)

        if response.status_code != 200:
            print(response.status_code)
            print(response.text)
            response.raise_for_status()

        data = response.json()

        results = data["data"]
        all_results.extend(results)

    return all_results



def get_objects(token, objectid_list, batch_size=1000):
    url = f"{BOOM_URL}/queries/find"

    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }

    all_results = []

    for i in range(0, len(objectid_list), batch_size):
        batch_ids = objectid_list[i:i + batch_size]

        payload = {
            "catalog_name": "LSST_alerts_aux",
            "filter": { "_id": {"$in": batch_ids} },
            "projection": {
                "objectId": 1,
                "prv_candidates.jd": 1,
                "prv_candidates.magpsf": 1,
                "prv_candidates.sigmapsf": 1,          
                "prv_candidates.band": 1,
                "prv_candidates.snr": 1
            }
        }

        print(f"Querying object batch {i // batch_size + 1} ({len(batch_ids)} ids)")

        response = requests.post(url, json=payload, headers=headers)

        if response.status_code != 200:
            print(response.status_code)
            print(response.text)
            response.raise_for_status()

        data = response.json()

        results = data["data"]
        all_results.extend(results)

    return all_results


# MAIN
# =========================
# INSPECT: NED + SIMBAD + ALADIN HTML
def _require_inspect_deps():
    try:
        from astropy.coordinates import SkyCoord
        from astroquery.simbad import Simbad
        return SkyCoord, Simbad
    except ImportError:
        raise ImportError(
            "INSPECT requires astroquery. Run:\n"
            "  pip install astroquery"
        )

def _inspect_query_ned(ra, dec, radius_arcsec, timeout=20):
    import requests as _req
    rows = []
    try:
        from astroquery.ipac.ned import Ned
        from astropy.coordinates import SkyCoord
        import astropy.units as _u
        coord = SkyCoord(ra=ra * _u.deg, dec=dec * _u.deg)
        result = Ned.query_region(coord, radius=radius_arcsec * _u.arcsec)
        if result is not None and len(result) > 0:
            cols = result.colnames
            sep_col = next((c for c in ("Separation", "Distance (arcmin)", "AngDist", "ang_dist") if c in cols), None)
            if sep_col is None:
                print(f"    [NED] unknown sep column: {cols}")
            else:
                sep_unit = 60.0 if "arcmin" in sep_col.lower() else 1.0
                name_col = next((c for c in ("Object Name", "Name") if c in cols), cols[0])
                type_col = next((c for c in ("Type", "Object Type") if c in cols), "")
                z_col    = next((c for c in ("Redshift", "Velocity") if c in cols), "")
                for row in result:
                    z_raw = row[z_col] if z_col else None
                    try:
                        z = float(z_raw) if z_raw is not None and str(z_raw) not in ("--","nan","","N/A") else None
                    except Exception:
                        z = None
                    rows.append({
                        "name":       str(row[name_col]),
                        "ra":         float(row["RA"]),
                        "dec":        float(row["DEC"]),
                        "type":       str(row[type_col]) if type_col else "?",
                        "redshift":   z,
                        "sep_arcsec": float(row[sep_col]) * sep_unit,
                    })
                return rows
    except Exception as e:
        print(f"    [NED] astroquery failed: {e} — trying REST fallback")

    #REST fallback
    radius_arcmin = radius_arcsec / 60.0
    url = "https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch"
    params = {
        "search_type": "Near Position Search",
        "in_csys": "Equatorial", "in_equinox": "J2000.0",
        "lon": f"{ra}d", "lat": f"{dec}d",
        "radius": str(radius_arcmin),
        "out_csys": "Equatorial", "out_equinox": "J2000.0",
        "of": "ascii_bar", "zv_breaker": "30000.0",
        "list_limit": "100", "img_stamp": "NO",
    }
    try:
        resp = _req.get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        for line in resp.text.splitlines():
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("No."):
                continue
            parts = line.split("|")
            if len(parts) < 13:
                continue
            try:
                z_raw = parts[6].strip()
                z = float(z_raw) if z_raw and z_raw not in ("N/A","Unknown","") else None
                sep_raw = parts[12].strip()
                sep_arcsec = float(sep_raw) * 60.0 if sep_raw else None
                if sep_arcsec is None:
                    continue
                rows.append({
                    "name":       parts[1].strip(),
                    "ra":         float(parts[2].strip()) if parts[2].strip() else ra,
                    "dec":        float(parts[3].strip()) if parts[3].strip() else dec,
                    "type":       parts[4].strip(),
                    "redshift":   z,
                    "sep_arcsec": sep_arcsec,
                })
            except (ValueError, IndexError):
                continue
        return rows
    except Exception as e:
        print(f"    [NED REST] also failed: {e}")
        return []


def _inspect_query_simbad(SkyCoord, ra, dec, radius_arcsec, timeout=20):
    import requests as _req
    rows = []

    radius_deg = round(radius_arcsec / 3600.0, 8)
    adql = (
        f"SELECT main_id, ra, dec, otype_txt, rvz_redshift "
        f"FROM basic "
        f"WHERE CONTAINS(POINT(ra,dec),CIRCLE({ra},{dec},{radius_deg}))=1 "
        f"LIMIT 100"
    )
    tap_url = "https://simbad.u-strasbg.fr/simbad/sim-tap/sync"
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": adql}
    try:
        resp = _req.get(tap_url, params=params, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
        col_names = [c["name"] for c in data.get("metadata", [])]
        for row_vals in data.get("data", []):
            row = dict(zip(col_names, row_vals))
            z_raw = row.get("rvz_redshift")
            try:
                z = float(z_raw) if z_raw is not None and str(z_raw) not in ("--","nan","") else None
            except Exception:
                z = None
            rows.append({
                "name":  str(row.get("main_id", "?")),
                "ra":    row.get("ra"),
                "dec":   row.get("dec"),
                "otype": str(row.get("otype_txt", "?")),
                "z":     z,
            })
        if rows:
            return rows
    except Exception as e:
        print(f"    [SIMBAD TAP] failed: {e} — trying astroquery fallback")

    #astroquery fallback. 
    try:
        import warnings, astropy.units as _u
        from astroquery.simbad import Simbad as _Simbad
        sim = _Simbad()
        sim.ROW_LIMIT = 100
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for field in ("otype", "rvz_redshift", "rvz_type", "rvz_qual", "z_value", "rvz_error"):
                try:
                    sim.add_votable_fields(field)
                except Exception:
                    pass
        coord = SkyCoord(ra=ra * _u.deg, dec=dec * _u.deg)
        result = sim.query_region(coord, radius=radius_arcsec * _u.arcsec)
        if result is None or len(result) == 0:
            return []
        col_upper = {c.upper(): c for c in result.colnames}
        def gcol(*names):
            for n in names:
                if n.upper() in col_upper:
                    return col_upper[n.upper()]
            return None
        name_col  = gcol("MAIN_ID", "main_id", "NAME")
        otype_col = gcol("OTYPE", "otype_txt")
        z_col     = gcol("RVZ_REDSHIFT", "Z_VALUE", "rvz_redshift", "z_value")
        ra_col    = gcol("RA_d", "RA", "ra_d")
        dec_col   = gcol("DEC_d", "DEC", "dec_d")
        if name_col is None:
            print(f"    [SIMBAD] no name column in: {result.colnames}")
            return []
        for row in result:
            z_raw = row[z_col] if z_col else None
            try:
                z = float(z_raw) if z_raw is not None and str(z_raw) not in ("--","nan","","N/A") else None
            except Exception:
                z = None
            rows.append({
                "name":  str(row[name_col]),
                "ra":    float(row[ra_col])  if ra_col  else None,
                "dec":   float(row[dec_col]) if dec_col else None,
                "otype": str(row[otype_col]) if otype_col else "?",
                "z":     z,
            })
        return rows
    except Exception as e:
        print(f"    [SIMBAD] astroquery fallback also failed: {e}")
        return []

def _inspect_validate_ned(boom_z, aq_ned_rows, sep_threshold=10.0, z_tol=0.05):
    empty = {
        "status": "no_aq_match", "aq_z": None, "delta_z": None,
        "matched_name": None, "matched_ra": None, "matched_dec": None,
        "matched_type": None, "matched_sep_arcsec": None,
    }
    if not aq_ned_rows:
        return empty

    close = [r for r in aq_ned_rows if r.get("sep_arcsec") is not None and r["sep_arcsec"] < sep_threshold]
    if not close:
        return empty

    def pack(status, row, dz=None):
        return {
            "status": status,
            "aq_z": row.get("redshift") if row else None,
            "delta_z": round(dz, 5) if dz is not None else None,
            "matched_name": row.get("name") if row else None,
            "matched_ra": row.get("ra") if row else None,
            "matched_dec": row.get("dec") if row else None,
            "matched_type": row.get("type") if row else None,
            "matched_sep_arcsec": row.get("sep_arcsec") if row else None,
        }

    with_z = [r for r in close if r.get("redshift") is not None]
    if not with_z:
        if boom_z is None:
            return {**empty, "status": "boom_missing_z"}
        return empty

    if boom_z is None:
        best = min(with_z, key=lambda r: r["sep_arcsec"])
        return pack("boom_missing_z", best)

    z_candidates = [r for r in with_z if abs(r["redshift"] - boom_z) < z_tol]
    if z_candidates:
        #prefer redshift agreement first, then angular separation
        best = min(z_candidates, key=lambda r: (abs(r["redshift"] - boom_z), r["sep_arcsec"]))
    else:
        #if no redshift agreement, expose the nearest redshift-bearing object
        best = min(with_z, key=lambda r: r["sep_arcsec"])

    dz = abs(boom_z - best["redshift"])
    return pack("agree" if dz < 0.01 else "disagree", best, dz)

def _inspect_dm_from_z(z):
    if z is None or z <= 0:
        return None
    try:
        dl_pc = COSMO.luminosity_distance(z).to(u.pc).value
        return 5.0 * (math.log10(dl_pc) - 1.0)
    except Exception:
        return None


def _kv(label, value, highlight=False):
    cls = ' class="highlight"' if highlight else ''
    return '<div class="kv"' + cls + '><span class="kv-label">' + label + '</span><span class="kv-value">' + str(value) + '</span></div>'


def _boom_ned_html(cand):
    if not cand.get("has_nearby_host"):
        return '<div class="no-match">No BOOM NED host within 10&Prime;</div>'
    parts = [
        _kv("z",         f'{cand["host_z"]:.5f}'                 if cand["host_z"]          else "&mdash;"),
        _kv("sep",       f'{cand["host_sep_arcsec"]:.2f}&Prime;'  if cand["host_sep_arcsec"] else "&mdash;"),
        _kv("sep kpc",   f'{cand["host_sep_kpc"]:.1f}'            if cand["host_sep_kpc"]    else "&mdash;"),
        _kv("z_tech",    str(cand.get("z_tech")  or "&mdash;")),
        _kv("z_qual",    str(cand.get("z_qual")  or "&mdash;")),
        _kv("DM",        f'{cand["dm"]:.3f}'                      if cand["dm"]              else "&mdash;"),
        _kv("DM source", str(cand.get("dm_source") or "&mdash;")),
        _kv("M (obs)",   f'{cand["abs_mag"]:.2f}'                 if cand["abs_mag"]         else "&mdash;"),
        _kv("M (corr)",  f'{cand["abs_mag_corr"]:.2f}'            if cand["abs_mag_corr"]    else "&mdash;", highlight=True),
        _kv("E(B-V)",    f'{cand["ebv"]:.4f}'                     if cand["ebv"]             else "&mdash;"),
        _kv("A_band",    f'{cand["a_band"]:.4f}'                  if cand["a_band"]          else "&mdash;"),
        _kv("low-z unc.", "&#9888; YES" if cand.get("lowz_uncertain") else "no"),
    ]
    return "\n".join(parts)


def _aq_ned_html(rows, radius):
    if not rows:
        return '<div class="no-match">No NED objects within ' + str(radius) + '&Prime;</div>'
    html = []
    for r in sorted(rows, key=lambda x: x["sep_arcsec"])[:10]:
        dm = _inspect_dm_from_z(r["redshift"])
        z_span  = '<span class="xm-z">z=' + f'{r["redshift"]:.5f}' + '</span>' if r["redshift"] else '<span class="xm-z muted">no z</span>'
        dm_span = '<span class="xm-dm">DM=' + f'{dm:.2f}' + '</span>' if dm else ''
        html.append(
            '<div class="xm-row">'
            '<span class="xm-name">' + r["name"] + '</span>'
            '<span class="xm-type">' + r["type"] + '</span>'
            '<span class="xm-sep">' + f'{r["sep_arcsec"]:.1f}' + '&Prime;</span>'
            + z_span + dm_span + '</div>'
        )
    if len(rows) > 10:
        html.append('<div class="muted">&hellip; ' + str(len(rows) - 10) + ' more</div>')
    return "\n".join(html)

def _simbad_html(rows):
    if not rows:
        return '<div class="no-match">No SIMBAD objects in cone</div>'
    html = []
    for r in rows[:10]:
        z_span = '<span class="xm-z">z=' + f'{r["z"]:.5f}' + '</span>' if r["z"] else '<span class="xm-z muted">no z</span>'
        html.append(
            '<div class="xm-row">'
            '<span class="xm-name">' + r["name"] + '</span>'
            '<span class="xm-type">' + r["otype"] + '</span>'
            + z_span + '</div>'
        )
    if len(rows) > 10:
        html.append('<div class="muted">&hellip; ' + str(len(rows) - 10) + ' more</div>')
    return "\n".join(html)

def _validation_html(val, boom_z):
    status = val["status"]
    color = {"agree": "#22c55e", "disagree": "#ef4444",
             "no_aq_match": "#94a3b8", "boom_missing_z": "#f97316"}.get(status, "#94a3b8")
    bz   = f'{boom_z:.5f}'         if boom_z                    else "&mdash;"
    aqz  = f'{val["aq_z"]:.5f}'    if val["aq_z"]               else "&mdash;"
    dz   = f'{val["delta_z"]:.5f}' if val["delta_z"] is not None else "&mdash;"
    name = str(val.get("matched_name") or "&mdash;")
    parts = [
        '<div class="val-status" style="color:' + color + '">&#9679; ' + status.replace("_", " ").upper() + '</div>',
        _kv("BOOM z",        bz),
        _kv("AQ NED z",      aqz),
        _kv("&Delta;z",      dz,   highlight=(val.get("delta_z") or 0) > 0.01),
        _kv("matched object", name),
    ]
    return "\n".join(parts)

def _html_escape(value):
    import html as _html
    if value is None:
        return "&mdash;"
    return _html.escape(str(value), quote=True)

def _external_links_html(ra, dec, radius_arcsec, fov_deg):
    import urllib.parse as _url
    coord = f"{ra} {dec}"
    radius_arcmin = radius_arcsec / 60.0
    aladin_url = (
        "https://aladin.cds.unistra.fr/AladinLite/?target="
        + _url.quote(coord)
        + f"&fov={fov_deg}&survey=P%2FDSS2%2Fred"
    )
    simbad_url = (
        "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord="
        + _url.quote_plus(coord)
        + f"&Radius={radius_arcsec}&Radius.unit=arcsec"
    )
    ned_url = (
        "https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?"
        "search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0"
        + f"&lon={ra}d&lat={dec}d&radius={radius_arcmin}"
        "&out_csys=Equatorial&out_equinox=J2000.0&of=pre_text&list_limit=100&img_stamp=YES"
    )
    legacy_url = f"https://www.legacysurvey.org/viewer?ra={ra}&dec={dec}&zoom=16&layer=ls-dr10"
    return (
        '<div class="inspect-tools">'
        '<a class="link-btn" target="_blank" rel="noopener" href="' + aladin_url + '">Open Aladin/CDS</a>'
        '<a class="link-btn" target="_blank" rel="noopener" href="' + simbad_url + '">Open SIMBAD cone</a>'
        '<a class="link-btn" target="_blank" rel="noopener" href="' + ned_url + '">Open NED cone</a>'
        '<a class="link-btn" target="_blank" rel="noopener" href="' + legacy_url + '">Open Legacy Survey</a>'
        '</div>'
    )


def _safe_float(x):
    try:
        if x is None:
            return None
        if str(x).strip() in ("", "--", "nan", "None"):
            return None
        return float(x)
    except Exception:
        return None

def _angular_sep_arcsec(ra1, dec1, ra2, dec2):
    """Small, dependency-free angular separation in arcsec."""
    ra1 = _safe_float(ra1); dec1 = _safe_float(dec1)
    ra2 = _safe_float(ra2); dec2 = _safe_float(dec2)
    if None in (ra1, dec1, ra2, dec2):
        return None
    r1 = math.radians(ra1); d1 = math.radians(dec1)
    r2 = math.radians(ra2); d2 = math.radians(dec2)
    cosang = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    cosang = max(-1.0, min(1.0, cosang))
    return math.degrees(math.acos(cosang)) * 3600.0

def _fmt_arcsec(x):
    x = _safe_float(x)
    if x is None:
        return "unknown separation"
    if x < 0.05:
        return "<0.05\""
    return f'{x:.2f}\"'

def _catalog_warning_match(item, keywords, star_like=False):
    c = item.get("cand", {})
    cand_ra = c.get("ra")
    cand_dec = c.get("dec")
    matches = []

    for r in item.get("aq_ned", []):
        name = str(r.get("name", "unknown"))
        typ = str(r.get("type", ""))
        blob = (name + " " + typ).lower()
        is_match = (typ.strip() == "*") if star_like else any(k in blob for k in keywords)
        if is_match:
            sep = _safe_float(r.get("sep_arcsec"))
            if sep is None:
                sep = _angular_sep_arcsec(cand_ra, cand_dec, r.get("ra"), r.get("dec"))
            matches.append({"name": name, "type": typ or "NED", "sep": sep, "source": "NED"})

    for r in item.get("simbad", []):
        name = str(r.get("name", "unknown"))
        typ = str(r.get("otype", ""))
        blob = (name + " " + typ).lower()
        is_match = (typ.strip() == "*") if star_like else any(k in blob for k in keywords)
        if is_match:
            sep = _angular_sep_arcsec(cand_ra, cand_dec, r.get("ra"), r.get("dec"))
            matches.append({"name": name, "type": typ or "SIMBAD", "sep": sep, "source": "SIMBAD"})

    if not matches:
        return None
    #prefer the closest object if separations exist; otherwise first returned.
    return min(matches, key=lambda m: m["sep"] if m["sep"] is not None else 1e99)


def _contamination_warnings(item, radius_arcsec):
    c = item["cand"]
    warnings = []
    host_radius = MAX_HOST_ARCSEC

    #host-offset warnings from BOOM/NED host values
    sep = _safe_float(c.get("host_sep_arcsec"))
    sep_kpc = _safe_float(c.get("host_sep_kpc"))
    if sep is not None and sep < 0.75:
        if sep_kpc is not None:
            warnings.append(
                f'Selected BOOM/NED host is within the {host_radius:.0f}\" host-search radius '
                f'at {_fmt_arcsec(sep)} from the source ({sep_kpc:.1f} kpc projected) — '
                'very small offset; check nuclear/AGN contamination.'
            )
        else:
            warnings.append(
                f'Selected BOOM/NED host is within the {host_radius:.0f}\" host-search radius '
                f'at {_fmt_arcsec(sep)} from the source — very small offset; check nuclear/AGN contamination.'
            )
    if sep is not None and sep > 5.0:
        if sep_kpc is not None:
            warnings.append(
                f'Selected BOOM/NED host is within the {host_radius:.0f}\" host-search radius '
                f'at {_fmt_arcsec(sep)} from the source ({sep_kpc:.1f} kpc projected) — '
                'large offset.'
            )
        else:
            warnings.append(
                f'Selected BOOM/NED host is within the {host_radius:.0f}\" host-search radius '
                f'at {_fmt_arcsec(sep)} from the source — large offset.'
            )
    if c.get("lowz_uncertain"):
        warnings.append("Low-redshift distance from z is uncertain because peculiar velocities may matter.")

    #nearby catalogue contaminants
    qso = _catalog_warning_match(item, ["qso", "quasar", "bl lac", "bllac"])
    if qso:
        warnings.append(
            f'{qso["source"]} {qso["type"]} object "{qso["name"]}" is within the '
            f'{radius_arcsec:.0f}\" inspection cone at {_fmt_arcsec(qso["sep"])} from the source — '
            'possible background contaminant.'
        )

    agn = _catalog_warning_match(item, ["seyfert", " sy", "sy1", "sy2", "agn"])
    if agn:
        warnings.append(
            f'{agn["source"]} {agn["type"]} object "{agn["name"]}" is within the '
            f'{radius_arcsec:.0f}\" inspection cone at {_fmt_arcsec(agn["sep"])} from the source — '
            'possible AGN contamination for a kilonova search.'
        )

    sn = _catalog_warning_match(item, [" sn", "sn*", "supernova"])
    if sn:
        warnings.append(
            f'{sn["source"]} {sn["type"]} object "{sn["name"]}" is within the '
            f'{radius_arcsec:.0f}\" inspection cone at {_fmt_arcsec(sn["sep"])} from the source — '
            'known/previous supernova-like object; check history.'
        )

    star = _catalog_warning_match(item, [], star_like=True)
    if star:
        warnings.append(
            f'{star["source"]} star-like object "{star["name"]}" is within the '
            f'{radius_arcsec:.0f}\" inspection cone at {_fmt_arcsec(star["sep"])} from the source — '
            'check stellar contamination.'
        )
    deduped = []
    for w in warnings:
        if w not in deduped:
            deduped.append(w)
    return deduped

def _inspection_notes_html(item, radius_arcsec):
    warnings = _contamination_warnings(item, radius_arcsec)

    if warnings:
        warning_html = '<ul class="warning-list">' + "".join('<li>' + _html_escape(w) + '</li>' for w in warnings) + '</ul>'
        return (
            '<div class="verdict warn">'
            '<div class="verdict-title">Inspection warnings</div>'
            '<div class="verdict-msg">Catalogue/context warnings found in the inspection cone.</div>'
            + warning_html +
            '</div>'
        )

    return (
        '<div class="verdict ok">'
        '<div class="verdict-title">Inspection notes</div>'
        '<div class="verdict-msg">No obvious catalogue contamination warning from the cone-search summary.</div>'
        '</div>'
    )

def _build_inspect_card(idx, item, fov_deg, radius_arcsec):
    c      = item["cand"]
    val    = item["validation"]
    ra     = c.get("ra")  or 0.0
    dec    = c.get("dec") or 0.0
    oid    = c.get("objectId", "&mdash;")
    candid = str(c.get("_id", "&mdash;"))
    band   = c.get("band") or "&mdash;"
    magpsf = f'{c["magpsf"]:.3f}' if c.get("magpsf") else "&mdash;"
    snr    = f'{c["snr_psf"]:.1f}' if c.get("snr_psf") else "&mdash;"
    jd     = f'{c["jd"]:.4f}'      if c.get("jd")     else "&mdash;"

    boom_html    = _boom_ned_html(c)
    aq_html      = _aq_ned_html(item["aq_ned"], radius_arcsec)
    sim_html     = _simbad_html(item["simbad"])
    notes_html   = _inspection_notes_html(item, radius_arcsec)
    links_html   = _external_links_html(ra, dec, radius_arcsec, fov_deg)

    ra_s  = str(ra)
    dec_s = str(dec)
    fov_s = str(fov_deg)
    radius_deg_s = str(round(10.0 / 3600.0, 8))

    oid_js   = json.dumps(str(oid))
    band_js  = json.dumps(str(band))
    mag_js   = json.dumps(str(magpsf))

    host_ra = val.get("matched_ra")
    host_dec = val.get("matched_dec")
    host_name = val.get("matched_name")
    host_z = val.get("aq_z")
    host_sep = val.get("matched_sep_arcsec")
    has_host_marker = host_ra is not None and host_dec is not None
    host_ra_s = str(host_ra) if has_host_marker else "null"
    host_dec_s = str(host_dec) if has_host_marker else "null"
    host_name_js = json.dumps(str(host_name or "NED cone match"))
    host_z_js = json.dumps("—" if host_z is None else f"{host_z:.5f}")
    host_sep_js = json.dumps("—" if host_sep is None else f"{host_sep:.2f} arcsec")

    js = (
        "<script>\n"
        "(function() {\n"
        "  function showAladin" + str(idx) + "() {\n"
        "    var box = document.getElementById('aladin-" + str(idx) + "');\n"
        "    if (!box) return;\n"
        "    box.innerHTML = '';\n"
        "    if (!window.A || !A.aladin) {\n"
        "      box.textContent = 'Aladin library did not load. Check internet access to aladin.cds.unistra.fr.';\n"
        "      box.style.padding = '18px';\n"
        "      box.style.color = '#ef4444';\n"
        "      box.style.fontFamily = 'monospace';\n"
        "      return;\n"
        "    }\n"
        "    try {\n"
        "      var aladin = A.aladin('#aladin-" + str(idx) + "', {\n"
        "        survey: 'P/DSS2/red',\n"
        "        target: '" + ra_s + " " + dec_s + "',\n"
        "        fov: " + fov_s + ",\n"
        "        showReticle: true,\n"
        "        showZoomControl: true,\n"
        "        showFullscreenControl: true,\n"
        "        showLayersControl: true,\n"
        "        showGotoControl: true\n"
        "      });\n"
        "      try {\n"
        "        var candCat = A.catalog({name:'LSST candidate', color:'#ff0000', shape:'cross', sourceSize:24, onClick:'showPopup'});\n"
        "        candCat.addSources([A.source(" + ra_s + ", " + dec_s + ", {name:" + oid_js + ", magpsf:" + mag_js + ", band:" + band_js + "})]);\n"
        "        aladin.addCatalog(candCat);\n"
        "      } catch(markerErr) { console.warn('Candidate marker failed:', markerErr); }\n"
        "      try {\n"
        "        if (" + ("true" if has_host_marker else "false") + ") {\n"
        "          var hostCat = A.catalog({name:'NED cone match', color:'#22c55e', shape:'circle', sourceSize:18, onClick:'showPopup'});\n"
        "          hostCat.addSources([A.source(" + host_ra_s + ", " + host_dec_s + ", {name:" + host_name_js + ", z:" + host_z_js + ", sep:" + host_sep_js + "})]);\n"
        "          aladin.addCatalog(hostCat);\n"
        "        }\n"
        "      } catch(hostErr) { console.warn('Host marker failed:', hostErr); }\n"
        "      try {\n"
        "        if (A.graphicOverlay && A.circle && aladin.addOverlay) {\n"
        "          var overlay = A.graphicOverlay({name:'10 arcsec host-search radius', color:'#facc15', lineWidth:4});\n"
        "          aladin.addOverlay(overlay);\n"
        "          if (overlay.addFootprints) {\n"
        "            overlay.addFootprints([A.circle(" + ra_s + ", " + dec_s + ", " + radius_deg_s + ", {color:'#facc15', lineWidth:4})]);\n"
        "          } else if (overlay.add) {\n"
        "            overlay.add(A.circle(" + ra_s + ", " + dec_s + ", " + radius_deg_s + ", {color:'#facc15', lineWidth:4}));\n"
        "          }\n"
        "        }\n"
        "      } catch(circleErr) { console.warn('Search-radius circle failed:', circleErr); }\n"
        "    } catch(e) {\n"
        "      console.error('Aladin card " + str(idx) + " failed:', e);\n"
        "      box.textContent = 'Aladin failed: ' + (e && e.message ? e.message : e);\n"
        "      box.style.padding = '18px';\n"
        "      box.style.color = '#ef4444';\n"
        "      box.style.fontFamily = 'monospace';\n"
        "    }\n"
        "  }\n"
        "  window.addEventListener('load', function() {\n"
        "    if (window.A && A.init && typeof A.init.then === 'function') {\n"
        "      A.init.then(showAladin" + str(idx) + ").catch(function(e) {\n"
        "        console.error('A.init failed:', e);\n"
        "        showAladin" + str(idx) + "();\n"
        "      });\n"
        "    } else {\n"
        "      setTimeout(showAladin" + str(idx) + ", 1000);\n"
        "    }\n"
        "  });\n"
        "})();\n"
        "</script>\n"
    )

    return (
        '\n<div class="card" id="card-' + str(idx) + '">\n'
        '  <div class="card-header">\n'
        '    <div class="card-title">\n'
        '      <span class="obj-id">' + _html_escape(oid) + '</span>\n'
        '      <span class="candid">candid: ' + _html_escape(candid) + '</span>\n'
        '    </div>\n'
        '    <div class="card-meta">\n'
        '      <span class="tag band">band: ' + _html_escape(band) + '</span>\n'
        '      <span class="tag mag">m = ' + _html_escape(magpsf) + '</span>\n'
        '      <span class="tag snr">SNR = ' + _html_escape(snr) + '</span>\n'
        '      <span class="tag jd">JD ' + _html_escape(jd) + '</span>\n'
        '    </div>\n'
        '  </div>\n'
        '  <div class="card-body">\n'
        '    <div class="aladin-wrap">\n'
        '      <div id="aladin-' + str(idx) + '" style="width:100%;height:420px;background:#05070a;"></div>\n'
        '      <div class="marker-legend"><span class="dot-red"></span>candidate <span class="dot-green"></span>NED cone match <span class="dot-yellow"></span>10&Prime;&Prime; search radius</div>\n'
        + notes_html + '\n'
        + links_html + '\n'
        '    </div>\n'
        '    <div class="xmatch-panels">\n'
        '      <div class="panel"><div class="panel-title">BOOM NED host</div>' + boom_html + '</div>\n'
        '      <div class="panel"><div class="panel-title">NED cone results <span class="radius-note">(' + str(radius_arcsec) + '&Prime;)</span></div>' + aq_html + '</div>\n'
        '      <div class="panel"><div class="panel-title">SIMBAD <span class="radius-note">(' + str(radius_arcsec) + '&Prime;)</span></div>' + sim_html + '</div>\n'
        '    </div>\n'
        '  </div>\n'
        '</div>\n'
        + js
    )

def _build_inspect_html(enriched, radius_arcsec, fov_arcmin):
    fov_deg = round(fov_arcmin / 60.0, 6)
    cards = [_build_inspect_card(idx, item, fov_deg, radius_arcsec) for idx, item in enumerate(enriched)]
    n         = len(enriched)
    with_host = sum(1 for d in enriched if d["cand"].get("has_nearby_host"))
    css = (
        "*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}\n"
        ":root{--bg:#0a0c10;--bg2:#10141c;--bg3:#161b26;--border:#1e2535;--accent:#38bdf8;"
        "--accent2:#f97316;--green:#22c55e;--red:#ef4444;--muted:#4b5563;--text:#e2e8f0;"
        "--text-dim:#94a3b8;--mono:'Space Mono',monospace;--sans:'Inter',sans-serif}\n"
        "body{background:var(--bg);color:var(--text);font-family:var(--sans);font-size:13px;line-height:1.5}\n"
        "header{background:var(--bg2);border-bottom:1px solid var(--border);padding:20px 32px;display:flex;align-items:center;gap:20px}\n"
        "header h1{font-family:var(--mono);font-size:16px;color:var(--accent);letter-spacing:.08em;text-transform:uppercase}\n"
        "header p{font-size:12px;color:var(--text-dim);font-family:var(--mono)}\n"
        ".dot{width:8px;height:8px;background:var(--red);border-radius:50%;animation:pulse 1.8s ease-in-out infinite;flex-shrink:0}\n"
        "@keyframes pulse{0%,100%{opacity:1;transform:scale(1)}50%{opacity:.4;transform:scale(1.4)}}\n"
        ".summary-bar{background:var(--bg3);border-bottom:1px solid var(--border);padding:12px 32px;display:flex;gap:24px;flex-wrap:wrap;font-family:var(--mono);font-size:12px}\n"
        ".sum-item{color:var(--text-dim)}.sum-item b{color:var(--text)}.sum-item.warn b{color:var(--accent2)}\n"
        ".cards-container{padding:24px 32px;display:flex;flex-direction:column;gap:32px}\n"
        ".card{background:var(--bg2);border:1px solid var(--border);border-radius:8px;overflow:hidden}\n"
        ".card-header{background:var(--bg3);border-bottom:1px solid var(--border);padding:14px 20px;display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:10px}\n"
        ".card-title{display:flex;align-items:baseline;gap:12px}\n"
        ".obj-id{font-family:var(--mono);font-size:15px;font-weight:700;color:var(--accent)}\n"
        ".candid{font-family:var(--mono);font-size:11px;color:var(--muted)}\n"
        ".card-meta{display:flex;gap:8px;flex-wrap:wrap}\n"
        ".tag{font-family:var(--mono);font-size:11px;padding:2px 8px;border-radius:4px;border:1px solid var(--border);color:var(--text-dim)}\n"
        ".tag.band{border-color:#38bdf840;color:var(--accent)}.tag.mag{border-color:#f9731640;color:var(--accent2)}\n"
        ".card-body{display:grid;grid-template-columns:minmax(420px,1fr) minmax(420px,1fr)}\n"
        ".aladin-wrap{border-right:1px solid var(--border);min-height:420px;background:#05070a}\n"
        ".aladin-wrap .aladin-container,.aladin-wrap canvas{max-width:100%}\n"
        ".marker-legend{font-family:var(--mono);font-size:11px;color:var(--text-dim);padding:8px 12px;border-top:1px solid var(--border);display:flex;gap:10px;align-items:center;flex-wrap:wrap}\n"
        ".dot-red,.dot-green,.dot-yellow{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:3px}.dot-red{background:#ff0000}.dot-green{background:#22c55e}.dot-yellow{background:#facc15}\n"
        ".inspect-tools{padding:10px 12px 14px;display:flex;gap:8px;flex-wrap:wrap;border-top:1px solid var(--border)}\n"
        ".link-btn{font-family:var(--mono);font-size:11px;text-decoration:none;color:var(--accent);border:1px solid #38bdf840;border-radius:4px;padding:5px 8px;background:#0f172a}\n"
        ".link-btn:hover{border-color:var(--accent);background:#111c2d}\n"
        ".verdict{margin:12px;border:1px solid var(--border);border-left-width:4px;border-radius:6px;padding:10px 12px;background:#0d111a}\n"
        ".verdict.ok{border-left-color:var(--green)}.verdict.warn{border-left-color:var(--accent2)}.verdict.bad{border-left-color:var(--red)}\n"
        ".verdict-title{font-family:var(--mono);font-weight:700;font-size:12px;margin-bottom:4px}.verdict.ok .verdict-title{color:var(--green)}.verdict.warn .verdict-title{color:var(--accent2)}.verdict.bad .verdict-title{color:var(--red)}\n"
        ".verdict-msg{font-size:12px;color:var(--text-dim);margin-bottom:6px}.warning-list{padding-left:18px;color:#fbbf24;font-size:12px}.warning-list li{margin:2px 0}.no-warning{font-size:12px;color:var(--text-dim)}\n"
        "@media(max-width:1000px){.card-body{grid-template-columns:1fr}.aladin-wrap{border-right:none;border-bottom:1px solid var(--border)}.xmatch-panels{grid-template-columns:1fr}}\n"
        ".xmatch-panels{display:grid;grid-template-columns:1fr 1fr;grid-template-rows:auto auto}\n"
        ".panel{padding:14px 16px;border-bottom:1px solid var(--border);border-right:1px solid var(--border)}\n"
        ".panel:nth-child(2n){border-right:none}.panel:nth-last-child(-n+2){border-bottom:none}\n"
        ".panel-title{font-family:var(--mono);font-size:11px;text-transform:uppercase;letter-spacing:.1em;color:var(--text-dim);margin-bottom:10px;padding-bottom:6px;border-bottom:1px solid var(--border)}\n"
        ".radius-note{color:var(--muted);font-size:10px}\n"
        ".kv{display:flex;justify-content:space-between;padding:2px 0;border-bottom:1px solid #1e253530}\n"
        ".kv:last-child{border-bottom:none}.kv-label{color:var(--text-dim);font-size:12px}.kv-value{font-family:var(--mono);font-size:12px;color:var(--text)}.kv.highlight .kv-value{color:var(--accent);font-weight:700}\n"
        ".xm-row{display:flex;align-items:baseline;gap:6px;padding:3px 0;border-bottom:1px solid #1e253530;flex-wrap:wrap;font-size:12px}\n"
        ".xm-row:last-child{border-bottom:none}\n"
        ".xm-name{font-family:var(--mono);color:var(--text);font-size:11px;flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}\n"
        ".xm-type{color:var(--accent2);font-size:11px;flex-shrink:0}.xm-sep{color:var(--text-dim);font-size:11px;flex-shrink:0}\n"
        ".xm-z{font-family:var(--mono);color:var(--accent);font-size:11px;flex-shrink:0}.xm-z.muted{color:var(--muted)}\n"
        ".xm-dm{font-family:var(--mono);color:var(--green);font-size:10px;flex-shrink:0}\n"
        ".no-match{color:var(--muted);font-style:italic;font-size:12px;padding:4px 0}\n"
        ".val-status{font-family:var(--mono);font-size:13px;font-weight:700;margin-bottom:8px}\n"
        ".muted{color:var(--muted);font-size:11px;margin-top:4px}\n"
    )
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\"/>\n"
        "<title>NED Crossmatched</title>\n"
        "<link href=\"https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=Inter:wght@300;400;500;600&display=swap\" rel=\"stylesheet\"/>\n"
        "<script src=\"https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js\" charset=\"utf-8\"></script>\n"
        "<style>\n" + css + "</style>\n</head>\n<body>\n"
        "<header>\n  <div class=\"dot\"></div>\n  <div>\n"
        "    <h1>NED Crossmatched</h1>\n"
        "    <p>NED + SIMBAD &middot; Aladin Lite v3 &middot; " + str(n) + " candidates &middot; cone " + str(radius_arcsec) + "&Prime;</p>\n"
        "  </div>\n</header>\n"
        "<div class=\"summary-bar\">\n"
        "  <span class=\"sum-item\">BOOM candidates with host crossmatched: <b>" + str(with_host) + "</b></span>\n"
        "</div>\n<div class=\"cards-container\">\n"
        + "".join(cards)
        + "\n</div>\n</body></html>"
    )


def run_inspect(unique, output_dir_flagged, savetime, range_string,
                radius_arcsec=60.0, fov_arcmin=3.0):
    """Re-query NED+SIMBAD and generate Aladin HTML for unique flagged objects.
    Called from main() only when INSPECT is True."""
    SkyCoord, Simbad = _require_inspect_deps()
    print(f"\n[inspect] Querying NED + SIMBAD for {len(unique)} unique candidates...")
    enriched = []
    for i, cand in enumerate(unique):
        ra  = cand.get("ra")
        dec = cand.get("dec")
        oid = cand.get("objectId", f"cand_{i}")
        print(f"  [{i+1}/{len(unique)}] {oid}  RA={ra}  Dec={dec}")
        if ra is None or dec is None:
            print("    WARNING: No RA/Dec - skipping")
            aq_ned, simbad = [], []
        else:
            aq_ned = _inspect_query_ned(ra, dec, radius_arcsec)
            print(f"    NED: {len(aq_ned)} objects")
            simbad = _inspect_query_simbad(SkyCoord, ra, dec, radius_arcsec)
            print(f"    SIMBAD: {len(simbad)} objects")
        val = _inspect_validate_ned(cand.get("host_z"), aq_ned, sep_threshold=10.0)
        print(f"    NED cone match for marker: {val.get('matched_name')}")
        enriched.append({"cand": cand, "aq_ned": aq_ned, "simbad": simbad, "validation": val})
    html = _build_inspect_html(enriched, radius_arcsec, fov_arcmin)
    out_path = os.path.join(output_dir_flagged, f"inspect_{savetime}_range_{range_string}.html")
    with open(out_path, "w") as f:
        f.write(html)
    print(f"[inspect] Saved to {out_path}\n")


def main(filter_file, output_dir_candidates_json, output_dir_objects_json, output_dir_summary, output_dir_flagged, start_time, end_time, save=False,
         inspect=False, inspect_radius=60.0, inspect_fov=3.0, skip_host=False):

    pipeline = load_pipeline(filter_file)

    print("Auth...")
    token, _ = get_token()

    start_dt = Time(start_time, scale='utc').to_datetime()
    end_dt = Time(end_time, scale='utc').to_datetime()

    all_alerts = []

    print("Running LSST filter in 1-day batches...\n")

    current = start_dt

    while current < end_dt:
        next_day = min(current + timedelta(days=1), end_dt)

        start_iso = current.isoformat()
        end_iso = next_day.isoformat()

        print(f"Querying: {start_iso} to {end_iso}")

        result = run_filter(token, pipeline, start_iso, end_iso)
        alerts = result.get("data", {}).get("results", [])

        print(f"  fetched {len(alerts)} alerts")

        all_alerts.extend(alerts)

        current = next_day

    print(f"\n                  Total alerts collected: {len(all_alerts)}")

    id_list = [a["_id"] for a in all_alerts]
    objectid_list = list(dict.fromkeys(a["objectId"] for a in all_alerts))

    print("===================================================================")
    print(f"                => Unique candidate IDs: {len(id_list)}")
    print(f"                => Unique object IDs: {len(objectid_list)}")
    print("===================================================================\n")

    candidate_jsons = get_candidates(token, id_list)
    ned_by_id = {}
    for a in all_alerts:
        ann = a.get("annotations", {}) or {}
        ned_by_id[a["_id"]] = ann.get("Cross-match_NED") or a.get("cross_matches", {}).get("NED") or []

    # HOST ASSOCIATION, ONE RUN PER OBJECT
    # =========================
    object_representatives = {}
    for candidate in candidate_jsons:
        cand = candidate.get("candidate", {})
        cand_ra = cand.get("ra")
        cand_dec = cand.get("dec")
        if cand_ra is None or cand_dec is None:
            continue

        obj_id = str(candidate.get("objectId", candidate.get("_id", "unknown")))
        snr = cand.get("snr_psf")
        mag = cand.get("magpsf")
        if snr is not None:
            score = float(snr)
        elif mag is not None:
            score = -float(mag)  #brighter is better if SNR is absent
        else:
            score = -1.0

        old_rep = object_representatives.get(obj_id)
        if old_rep is None or score > old_rep["score"]:
            object_representatives[obj_id] = {
                "ra": float(cand_ra),
                "dec": float(cand_dec),
                "score": score,
                "candid": candidate.get("_id"),
            }

    host_by_object = {}
    n_obj = len(object_representatives)
    savetime = Time.now().iso.replace(':', '-').replace(' ', '_')
    range_string = f"{start_time.replace(':','-')}_to_{end_time.replace(':','-')}"
    os.makedirs(output_dir_summary, exist_ok=True)
    incremental_csv = os.path.join(
        output_dir_summary,
        f"cands_incremental_{savetime}_range_{range_string}.csv"
    )
    if not skip_host:
        from host_association import associate_host
    if skip_host:
        print(f"\nSkipping host association (SKIP_HOST_ASSOCIATION): {len(object_representatives)} objects, host columns left empty.\n")
    else:
        print(f"\nRunning host association for {len(object_representatives)} unique objects, not {len(candidate_jsons)} candidate alerts.\n")
    for i, (obj_id, rep) in enumerate(object_representatives.items() if not skip_host else [], start=1):
        host_by_object[obj_id] = associate_host(
            transient_name=obj_id,
            ra_deg=rep["ra"],
            dec_deg=rep["dec"],
            obj_idx=i,
            obj_total=n_obj,
        )
        _partial_rows = []
        for _cand in candidate_jsons:
            _oid = str(_cand.get("objectId", _cand.get("_id", "unknown")))
            if _oid not in host_by_object:
                continue
            _hi = host_by_object[_oid]
            _c  = _cand.get("candidate", {})
            _partial_rows.append({k: _clean_for_output(_hi.get(k) if k in _hi else _cand.get(k) if k in _cand else _c.get(k)) for k, _ in CSV_COLUMNS})
        if _partial_rows:
            write_csv(incremental_csv, _partial_rows)

    #SUMMARY
    # =========================
    candidate_list = []

    for candidate in candidate_jsons:
        cand = candidate.get("candidate", {})
        magpsf  = cand.get("magpsf")
        cand_ra  = cand.get("ra")
        cand_dec = cand.get("dec")
        obj_id   = str(candidate.get("objectId", candidate.get("_id", "unknown")))
        host_info = host_by_object.get(obj_id, {})

        has_nearby_host = host_info.get("has_nearby_host", False)
        host_source      = host_info.get("host_source")
        host_z          = host_info.get("host_z")
        host_sep        = host_info.get("host_sep_arcsec")
        host_kpc        = host_info.get("host_sep_kpc")
        dm              = host_info.get("dm")
        dm_source       = host_info.get("dm_source")
        z_tech          = host_info.get("z_tech")
        z_unc           = host_info.get("z_unc")
        z_qual          = host_info.get("z_qual")

        ned  = ned_by_id.get(candidate.get("_id"), [])
        ebv  = None
        for m in ned:
            if m.get("ebv") is not None:
                ebv = m.get("ebv")
                break

        abs_mag = abs_mag_corr = a_band = None
        lowz_uncertain = False
        band = get_band(cand)

        if magpsf is not None and dm is not None:
            abs_mag = magpsf - dm

        a_band = galactic_extinction(ebv, band)
        if abs_mag is not None and a_band is not None:
            abs_mag_corr = abs_mag - a_band

        if dm_source is not None and "z_cosmo" in str(dm_source) and host_z is not None and host_z < LOWZ_UNCERTAIN:
            lowz_uncertain = True

        candidate_list.append({
            "_id":              candidate.get("_id"),
            "objectId":         candidate.get("objectId"),
            "ra":               cand_ra,
            "dec":              cand_dec,
            "jd":               cand.get("jd"),
            "band":             band,
            "magpsf":           magpsf,
            "snr_psf":          cand.get("snr_psf"),
            "has_nearby_host":  has_nearby_host,
            "host_source":       host_source,
            "host_z":           host_z,
            "z_unc":            z_unc,
            "z_tech":           z_tech,
            "z_qual":           z_qual,
            "host_sep_arcsec":  host_sep,
            "host_sep_kpc":     host_kpc,
            "ebv":              ebv,
            "a_band":           a_band,
            "dm":               dm,
            "dm_source":        dm_source,
            "abs_mag":          abs_mag,
            "abs_mag_corr":     abs_mag_corr,
            "lowz_uncertain":   lowz_uncertain,
            # GHOST
            "ghost_status":          host_info.get("ghost_status"),
            "ghost_ra":              host_info.get("ghost_ra"),
            "ghost_dec":             host_info.get("ghost_dec"),
            "ghost_z":               host_info.get("ghost_z"),
            "ghost_sep_arcsec":      host_info.get("ghost_sep_arcsec"),
            "ghost_objid":           host_info.get("ghost_objid"),
            "ghost_type":            host_info.get("ghost_type"),
            "ghost_error":           host_info.get("ghost_error"),
            # NED-LVS DLR
            "nedlvs_status":         host_info.get("nedlvs_status"),
            "nedlvs_ra":             host_info.get("nedlvs_ra"),
            "nedlvs_dec":            host_info.get("nedlvs_dec"),
            "nedlvs_dlr":            host_info.get("nedlvs_dlr"),
            "nedlvs_sep_arcsec":     host_info.get("nedlvs_sep_arcsec"),
            "nedlvs_objname":        host_info.get("nedlvs_objname"),
            "nedlvs_objtype":        host_info.get("nedlvs_objtype"),
            "nedlvs_z":              host_info.get("nedlvs_z"),
            "nedlvs_distmpc":        host_info.get("nedlvs_distmpc"),
            "nedlvs_n_candidates":   host_info.get("nedlvs_n_candidates"),
            "nedlvs_n_ranked":       host_info.get("nedlvs_n_ranked"),
            "nedlvs_diam_arcsec":    host_info.get("nedlvs_diam_arcsec"),
            "nedlvs_a_arcsec":       host_info.get("nedlvs_a_arcsec"),
            "nedlvs_b_arcsec":       host_info.get("nedlvs_b_arcsec"),
            "nedlvs_axis_ratio":     host_info.get("nedlvs_axis_ratio"),
            "nedlvs_pa_deg":         host_info.get("nedlvs_pa_deg"),
            # Agreement (all three pairwise separations)
            "host_agreement":            host_info.get("host_agreement"),
            "ghost_nedlvs_sep_arcsec":   host_info.get("ghost_nedlvs_sep_arcsec"),
            "ghost_legacy_sep_arcsec":   host_info.get("ghost_legacy_sep_arcsec"),
            "nedlvs_legacy_sep_arcsec":  host_info.get("nedlvs_legacy_sep_arcsec"),
            # Legacy TAP DLR
            "legacy_status":       host_info.get("legacy_status"),
            "legacy_ra":           host_info.get("legacy_ra"),
            "legacy_dec":          host_info.get("legacy_dec"),
            "legacy_dlr":          host_info.get("legacy_dlr"),
            "legacy_sep_arcsec":   host_info.get("legacy_sep_arcsec"),
            "legacy_type":         host_info.get("legacy_type"),
            "legacy_a_arcsec":     host_info.get("legacy_a_arcsec"),
            "legacy_b_arcsec":     host_info.get("legacy_b_arcsec"),
            "legacy_n_candidates": host_info.get("legacy_n_candidates"),
            "legacy_n_ranked":     host_info.get("legacy_n_ranked"),
        })

    output_csv = os.path.join(
        output_dir_summary,
        f"cands_{savetime}_range_{range_string}.csv"
    )

    write_csv(output_csv, candidate_list)

    if skip_host:
        print("\n                  SKIP_HOST_ASSOCIATION: no flagged / unique-flagged / inspect outputs.")
    else:
        #FLAGGED HOST-ASSOCIATED OUTPUT
        # =========================
        flagged = [c for c in candidate_list if c["has_nearby_host"]]

        print(f"\n                  Host-associated candidates with usable DM: {len(flagged)} / {len(candidate_list)}")

        os.makedirs(output_dir_flagged, exist_ok=True)

        flagged_csv = os.path.join(
            output_dir_flagged,
            f"flagged_{savetime}_range_{range_string}.csv"
        )

        write_csv(flagged_csv, flagged)

        flagged_json = os.path.join(
            output_dir_flagged,
            f"flagged_{savetime}_range_{range_string}.json"
        )
        with open(flagged_json, "w") as f:
            json.dump(flagged, f, indent=2, default=_json_default)

        print(f"Saved flagged list to {flagged_csv}")
        print(f"Saved flagged JSON to {flagged_json}")

        #UNIQUE FLAGGED OBJECTS (one row per objectId, at peak/brightest epoch)
        # =========================
        by_obj = {}
        for c in flagged:
            oid = c["objectId"]
            if oid not in by_obj or (c["magpsf"] is not None and c["magpsf"] < by_obj[oid]["magpsf"]):
                by_obj[oid] = c

        unique = list(by_obj.values())

        print(f"                  Unique flagged objects: {len(unique)}")
        for c in unique:
            am = f"{c['abs_mag']:.2f}" if c['abs_mag'] is not None else "NA"
            amc = f"{c['abs_mag_corr']:.2f}" if c['abs_mag_corr'] is not None else "NA"
            kpc = f"{c['host_sep_kpc']:.1f}" if c['host_sep_kpc'] is not None else "NA"
            zstr = f"{c['host_z']:.4f}" if c.get('host_z') is not None else "NA"
            sepstr = f"{c['host_sep_arcsec']:.2f}" if c.get('host_sep_arcsec') is not None else "NA"
            print(f"                    {c['objectId']}  source={c.get('host_source')}  z={zstr}  "
                  f"sep={sepstr}\" ({kpc} kpc)  "
                  f"M={am}  M_corr={amc}  ({c['dm_source']})")

        unique_csv = os.path.join(
            output_dir_flagged,
            f"flagged_unique_{savetime}_range_{range_string}.csv"
        )
        write_csv(unique_csv, unique)

        unique_json = os.path.join(
            output_dir_flagged,
            f"flagged_unique_{savetime}_range_{range_string}.json"
        )
        with open(unique_json, "w") as f:
            json.dump(unique, f, indent=2, default=_json_default)

        print(f"Saved unique flagged list to {unique_csv}")
        print(f"Saved unique flagged JSON to {unique_json}")

        if inspect:
            run_inspect(
                unique,
                output_dir_flagged,
                savetime,
                range_string,
                radius_arcsec=inspect_radius,
                fov_arcmin=inspect_fov,
            )

    #SAVE 
    if save:

        object_jsons = get_objects(token, objectid_list)
        print("")

        os.makedirs(output_dir_candidates_json, exist_ok=True)
        os.makedirs(output_dir_objects_json, exist_ok=True)

        output_candidates_json = os.path.join(
            output_dir_candidates_json,
            f"cands_{savetime}_range_{range_string}.json"
        )

        output_objects_json = os.path.join(
            output_dir_objects_json,
            f"objs_{savetime}_range_{range_string}.json"
        )

        with open(output_candidates_json, "w") as f:
            json.dump(candidate_jsons, f, indent=2, default=_json_default)
        print(f"Saved candidates to {output_candidates_json}")

        with open(output_objects_json, "w") as f:
            json.dump(object_jsons, f, indent=2, default=_json_default)
        print(f"Saved objects to {output_objects_json}")
    
    print(f"Saved summary to {output_csv}")



if __name__ == "__main__":
    import ghost_inputs as inputs

    if not os.path.isfile(inputs.FILTER_FILE):
        raise FileNotFoundError(f"FILTER_FILE not found: {inputs.FILTER_FILE}")
    if Time(inputs.END_TIME, scale="utc") <= Time(inputs.START_TIME, scale="utc"):
        raise ValueError("END_TIME must be after START_TIME")

    main(
        inputs.FILTER_FILE,
        os.path.join(inputs.OUTPUT_DIR, "candidate_jsons"),
        os.path.join(inputs.OUTPUT_DIR, "object_jsons"),
        os.path.join(inputs.OUTPUT_DIR, "candidate_quicklooks"),
        os.path.join(inputs.OUTPUT_DIR, "ned_flagged"),
        inputs.START_TIME,
        inputs.END_TIME,
        inputs.SAVE_JSON,
        inspect=inputs.INSPECT,
        inspect_radius=inputs.INSPECT_RADIUS_ARCSEC,
        inspect_fov=inputs.INSPECT_FOV_ARCMIN,
        skip_host=inputs.SKIP_HOST_ASSOCIATION,
    )
