import json
import os
import argparse
import requests
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from astropy.time import Time


# Color mapping
COLORS = {
    "u": "purple",
    "g": "green",
    "r": "red",
    "i": "orange",
    "z": "brown",
    "y": "gray"
}

BANDS = ["u", "g", "r", "i", "z", "y"]

# =========================
# CONFIG
# =========================

BOOM_URL = "https://api.kaboom.caltech.edu"
INCLUDE_ZTF = False
INCLUDE_FP = True

# =========================
# AUTH
# =========================
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

# =========================
# DATA LOADING
# =========================
def load_json(file_path):
    with open(file_path, "r") as f:
        return json.load(f)


def first_valid_coordinates(points):
    for point in points:
        ra = point.get("ra")
        dec = point.get("dec")
        if ra is not None and dec is not None:
            return ra, dec
    return None, None

# =========================
# ZTF CONE SEARCH
# =========================

def angular_sep(ra1, dec1, ra2, dec2):
    # small-angle approximation
    dra = (ra1 - ra2) * np.cos(np.radians(dec1))
    ddec = dec1 - dec2
    return np.sqrt(dra**2 + ddec**2)


def match_ztf_object(token, ra, dec, radius_arcsec=1.0, limit=100):
    url = f"{BOOM_URL}/queries/cone_search"

    radius_deg = radius_arcsec / 3600.0

    payload = {
        "catalog_name": "ZTF_alerts_aux",
        "projection": {
                "_id": 1,                       # ZTF ID of the object

                "prv_candidates.jd": 1,         # Photometry for the 5-sigma detections
                "prv_candidates.magpsf": 1,
                "prv_candidates.sigmapsf": 1,
                "prv_candidates.band": 1,
                "prv_candidates.snr_psf": 1,
                "prv_candidates.ra": 1,
                "prv_candidates.dec": 1,

                "prv_nondetections.jd": 1,     # Photometry for the non-detections (upper limits)
                "prv_nondetections.diffmaglim": 1,
                "prv_nondetections.band": 1,

        },
        "object_coordinates": {
            "additionalProperty": [ra, dec]
        },
        "radius": radius_deg,
        "unit": "Degrees",
        "limit": limit,
    }

    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }

    resp = requests.post(url, json=payload, headers=headers)

    if resp.status_code != 200:
        print(resp.status_code)
        print(resp.text)
        resp.raise_for_status()

    data = resp.json()["data"]["additionalProperty"]

    if len(data) == 0:
        print("No ZTF match found")
        return None, [], []

    # -----------------------------------
    # Find closest object
    # -----------------------------------
    best_obj = None
    best_dist = np.inf

    for obj in data:
        # get radec from first candidate
        candidates = obj.get("prv_candidates", [])
        if len(candidates) == 0:
            continue

        p0 = candidates[0]
        ra_ztf = p0.get("ra")
        dec_ztf = p0.get("dec")

        if ra_ztf is None or dec_ztf is None:
            continue

        print(f"Checking ZTF object {obj.get('_id')}: RA={ra_ztf:.4f}, Dec={dec_ztf:.4f} ...")

        d = angular_sep(ra, dec, ra_ztf, dec_ztf)

        print(f"  Angular separation: {d*3600:.2f} arcsec")

        if d < best_dist:
            best_dist = d
            best_obj = obj

    if best_obj is None:
        print("No valid ZTF object with coordinates")
        return None, [], []

    # -----------------------------------
    # Extract detections
    # -----------------------------------
    detections = []
    for p in best_obj.get("prv_candidates", []):
        if is_valid_point(p):
            detections.append(p)

    # -----------------------------------
    # Extract upper limits
    # -----------------------------------
    limits = []
    for p in best_obj.get("prv_nondetections", []):
        if is_valid_ulim(p):
            limits.append(p)

    print(f"Matched ZTF object: {best_obj.get('_id')} (dist ≈ {best_dist*3600:.2f} arcsec)")

    return best_obj, detections, limits

# =========================
# ALERT HIGHLIGHTING
# =========================
def get_alert_points_passed(alerts_data):
    """
    Build a mapping from objectId -> set of (band, jd) that passed the filter.
    """
    alert_points_passed = defaultdict(set)
    for alert in alerts_data:
        obj_id = alert["objectId"]
        band = alert["candidate"]["band"]
        jd = alert["candidate"]["jd"]
        
        if band is None or jd is None:
            continue

        alert_points_passed[obj_id].add((band, jd))
    return alert_points_passed

def is_highlight(band, jd, passed_points, tol=1e-6):
    for b, j in passed_points:
        if b == band and abs(jd - j) < tol:
            return True
    return False

# ====================================
# PHOTOMETRY QUALITY CHECK
# ====================================

def is_valid_point(p):
    bool_magpsf = "magpsf" in p and p["magpsf"] is not None
    bool_sigmapsf = "sigmapsf" in p and p["sigmapsf"] is not None
    bool_jd = "jd" in p and p["jd"] is not None
    bool_band = "band" in p and p["band"] is not None
    bool_snr = "snr_psf" in p and p["snr_psf"] is not None and p["snr_psf"] >= 2

    return (bool_magpsf and bool_sigmapsf and bool_jd and bool_band and bool_snr)

def is_valid_ulim(p):
    bool_jd = "jd" in p and p["jd"] is not None
    bool_lim = "diffmaglim" in p and p["diffmaglim"] is not None
    bool_band = "band" in p and p["band"] is not None

    return (bool_jd and bool_lim and bool_band)

# =========================
# LIGHT CURVE PLOTTING
# =========================
def plot_object_lc(token, obj, alert_points_passed, output_dir, iter):
    obj_id = obj.get("_id", obj.get("objectId"))

    print(f"\n\n\n      ========================== Object {obj_id} ({iter+1}) ==========================\n")

    lc_points = obj.get("prv_candidates", [])
    if INCLUDE_FP:
        lc_point_fp = obj.get("fp_hists", [])
    else:
        lc_point_fp = []

    obj_ra, obj_dec = first_valid_coordinates(lc_points)
    if obj_ra is not None and obj_dec is not None:
        print(f"Object coordinates: RA={obj_ra:.4f}, Dec={obj_dec:.4f}")
    else:
        print("Object coordinates unavailable; skipping ZTF crossmatch for this object")

    lc_points_ztf = []
    ztf_ulims = []
    if INCLUDE_ZTF and obj_ra is not None and obj_dec is not None:
        best_match, lc_points_ztf, ztf_ulims = match_ztf_object(token, ra=obj_ra, dec=obj_dec)

    if not lc_points and not lc_point_fp and not lc_points_ztf:
        print("No plottable photometry for this object")
        return 0

    # Organize points by band
    points_by_band = defaultdict(list)
    for p in lc_points:
        if p.get("band") is not None and p.get("jd") is not None and p.get("magpsf") is not None:
            points_by_band[p["band"]].append(p)
    
    points_by_band_fp = defaultdict(list)
    if INCLUDE_FP:
        for p in lc_point_fp:
            if is_valid_point(p):
                points_by_band_fp[p["band"]].append(p)
    
    points_by_band_ztf = defaultdict(list)
    if INCLUDE_ZTF:
        for p in lc_points_ztf:
            if is_valid_point(p):
                points_by_band_ztf[p["band"]].append(p)

    # Compute earliest JD across all bands
    all_plot_jds = (
        [p["jd"] for points in points_by_band.values() for p in points]
        + [p["jd"] for points in points_by_band_fp.values() for p in points]
        + [p["jd"] for points in points_by_band_ztf.values() for p in points]
    )
    if not all_plot_jds:
        print("No valid photometry survived quality cuts for this object")
        return 0

    jd_min = min(all_plot_jds)
    jd_min_iso = Time(jd_min, format='jd', scale='utc').to_datetime()

    # Setup subplots
    fig, axes = plt.subplots(3, 2, sharex=True, figsize=(12, 8))
    axes = axes.flatten()

    delta_list = [0 for band in BANDS]
    x_min, x_max = 1e30, 0

    highlighted_count = 0

    for i, band in enumerate(BANDS):
        ax = axes[i]
        band_points = points_by_band.get(band, [])
        band_points_fp = points_by_band_fp.get(band, [])
        band_points_ztf = points_by_band_ztf.get(band, [])
        band_color = COLORS.get(band, "black")

        ax.set_ylabel(f"{band}-band")
        if i >= 4:  # bottom row shows x-axis
            ax.set_xlabel("Age (days)")

        if not band_points and not band_points_fp and not band_points_ztf:
            continue
        
        if INCLUDE_FP:
            if not band_points_fp:
                print(f"No forced photometry points in band {band}")

        # ========================================
        #           Sort points by JD 
        # ========================================

        # LSST points
        band_points.sort(key=lambda x: x["jd"])
        x = [p["jd"] - jd_min for p in band_points]
        y = [p["magpsf"] for p in band_points]
        yerr = [p.get("sigmapsf", 0.0) for p in band_points]

        # Forced photometry points
        if INCLUDE_FP:
            band_points_fp.sort(key=lambda x: x["jd"])
            x_fp = [p["jd"] - jd_min for p in band_points_fp]
            y_fp = [p["magpsf"] for p in band_points_fp]
            yerr_fp = [p.get("sigmapsf", 0.0) for p in band_points_fp]
        else:
            x_fp, y_fp, yerr_fp = [], [], []

        # ZTF points
        if INCLUDE_ZTF:
            band_points_ztf.sort(key=lambda x: x["jd"])
            x_ztf = [p["jd"] - jd_min for p in band_points_ztf]
            y_ztf = [p["magpsf"] for p in band_points_ztf]
            yerr_ztf = [p.get("sigmapsf", 0.0) for p in band_points_ztf]
        else:
            x_ztf, y_ztf, yerr_ztf = [], [], []

        # Computing the maximum time range
        x_all = x + x_fp + x_ztf
        if x_all:
            x_min = min(x_min, min(x_all))
            x_max = max(x_max, max(x_all))

        # ========================================
        #      Highlighting relevant points
        # ========================================

        # Identify which points to highlight ( <=> the LSST points that passed the filter)
        highlights = [is_highlight(band, p["jd"], alert_points_passed[obj_id]) for p in band_points]

        # Updating the count of highlighted points
        highlighted_count += sum(highlights)

        # ========================================
        #      Plotting the light curves
        # ========================================

        # Plot forced photometry points (if available)
        if INCLUDE_FP and x_fp:
            if all(yerr_fp):
                ax.errorbar(x_fp, y_fp, yerr=yerr_fp, fmt='x', color=band_color, alpha=0.2, label='Forced photometry')
            else:
                ax.scatter(x_fp, y_fp, color=band_color, alpha=0.2, label='Forced photometry')

        # Plot ZTF points (if available)
        if INCLUDE_ZTF and x_ztf:
            if all(yerr_ztf):
                ax.errorbar(x_ztf, y_ztf, yerr=yerr_ztf, fmt='s', color='cyan', alpha=0.7, label='ZTF detections')
            else:
                ax.scatter(x_ztf, y_ztf, color='cyan', alpha=0.7, label='ZTF detections')

        # Plot normal LSST points
        x_norm = [xx for xx, h in zip(x, highlights) if not h]
        y_norm = [yy for yy, h in zip(y, highlights) if not h]
        yerr_norm = [ee for ee, h in zip(yerr, highlights) if not h]
        if x_norm:
            if all(yerr_norm):
                ax.errorbar(x_norm, y_norm, yerr=yerr_norm, fmt='o', color=band_color, alpha=0.7)             # label='Other points'
            else:
                ax.scatter(x_norm, y_norm, color=band_color, alpha=0.7)                                       # label='Other points'

        # Plot highlighted LSST points
        x_high = [xx for xx, h in zip(x, highlights) if h]
        y_high = [yy for yy, h in zip(y, highlights) if h]
        yerr_high = [ee for ee, h in zip(yerr, highlights) if h]

        delta = 0

        if x_high:
            ax.errorbar(x_high, y_high, yerr=yerr_high, fmt='o', color='black', markersize=8, markeredgecolor='yellow', markeredgewidth=1.5, label='Passed filter')
        
            # Shade the window that passed the cuts

            jd_start_rel = x_high[0]
            jd_stop_rel = x_high[-1]

            if len(x_high) >= 2:
                delta = (jd_stop_rel - jd_start_rel)
                ax.axvspan(jd_start_rel, jd_stop_rel, color="gray", alpha=0.5, label="Period of interest $dt_i$")

        ax.invert_yaxis()
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(loc='upper right')
        delta_list[i] = delta
    
    # Compute the overall period of interest ratio across all bands
    baseline = x_max - x_min
    object_delta = max(delta_list) / baseline if baseline > 0 else 0

    if object_delta > 0.1:
        print(f"Significant period of interest: Delta = {object_delta:.2f}")

    if object_delta == 0:
        object_delta = "single point"

    plt.suptitle(f"Light Curve for Object {obj_id} ({iter+1}), T0 = {jd_min_iso.strftime('%Y-%m-%d %H:%M:%S')} UTC, $\\Delta = {object_delta}$")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    output_path = os.path.join(output_dir, f"lc_{iter+1}_object_{obj_id}.png")
    plt.savefig(output_path)
    plt.close(fig)
    print(f"\nSaved to {output_path}")

    return highlighted_count


# =========================
# MAIN
# =========================
def main(objects_json_file, alerts_json_file, output_dir, include_ztf=None, include_fp=None, output_subdir=None, max_objects=None):
    global INCLUDE_ZTF, INCLUDE_FP

    if include_ztf is not None:
        INCLUDE_ZTF = include_ztf
    if include_fp is not None:
        INCLUDE_FP = include_fp

    os.makedirs(output_dir, exist_ok=True)
    token = None
    if INCLUDE_ZTF:
        token, expires_at = get_token()

    highlight_count = 0

    # Output subfolder setup
    if output_subdir is None:
        objs_file = os.path.basename(objects_json_file)
        stem = os.path.splitext(objs_file)[0]
        output_subdir = stem[5:] if stem.startswith("objs_") else stem
    output_dir = os.path.join(output_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)

    # Load local data
    objects_data = load_json(objects_json_file)
    alerts_data = load_json(alerts_json_file)

    # Alerts that passed the filter
    alert_points_passed = get_alert_points_passed(alerts_data)

    # Plotting loop
    if max_objects is not None:
        objects_data = objects_data[:max_objects]

    for i, obj in enumerate(objects_data):
        highlighted_count = plot_object_lc(token, obj, alert_points_passed, output_dir, i)
        highlight_count += highlighted_count

    print(f"Total highlighted points: {highlight_count}")
    return {
        "plots_dir": output_dir,
        "n_objects_plotted": len(objects_data),
        "n_highlighted_points": highlight_count,
        "include_ztf": INCLUDE_ZTF,
        "include_fp": INCLUDE_FP,
    }

"""
example command to run:
python3 plot_lc.py --objects_json path/to/objects_jsons/objects_2026-02-26T12-00-00_to_2026-03-08T12-00-00.json --alerts_json path/to/candidate_jsons/candidates_2026-02-26T12-00-00_to_2026-03-08T12-00-00.json --output_dir path/to/lc_plots
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot LSST light curves with highlights")
    parser.add_argument("--objects_json", required=True, help="Path to objects JSON file")
    parser.add_argument("--alerts_json", required=True, help="Path to alerts JSON file")
    parser.add_argument("--output_dir", default="./lc_plots/", help="Directory to save plots")
    parser.add_argument("--include_ztf", action="store_true", help="Include ZTF crossmatch photometry when available")
    parser.add_argument("--no_fp", action="store_true", help="Do not include LSST forced photometry")
    parser.add_argument("--output_subdir", default=None, help="Optional plot subfolder name")
    parser.add_argument("--max_objects", type=int, default=None, help="Optional maximum number of objects to plot")
    args = parser.parse_args()

    main(
        args.objects_json,
        args.alerts_json,
        args.output_dir,
        include_ztf=args.include_ztf,
        include_fp=not args.no_fp,
        output_subdir=args.output_subdir,
        max_objects=args.max_objects,
    )
