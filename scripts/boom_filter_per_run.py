import requests
from datetime import datetime, timedelta, timezone
from astropy.time import Time
import os
import json
import csv
from collections import Counter

"""
Example usage:

python3 path/to/boom_filter_per_run.py --output_dir path/to/boom_runs --filter_file path/to/filter.json --start_time 2026-02-26T12:00:00 --end_time 2026-03-08T12:00:00

This creates one run folder inside --output_dir. The run folder contains:

    cands.json
    objs.json
    cands_sum.csv
    objs_sum.csv
"""

# =========================
# CONFIG
# =========================

BOOM_URL = "https://api.kaboom.caltech.edu"

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
# LSST PIPELINE
# =========================

def load_pipeline(filter_file):
    with open(filter_file, "r") as f:
        pipeline = json.load(f)

    if not isinstance(pipeline, list):
        raise ValueError("Filter file must contain a JSON list (Mongo pipeline)")

    return pipeline


# =========================
# BOOM QUERIES
# =========================


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

    # Loop over batches
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
                "objectId": 1,                  # ID of the object

                "prv_candidates.jd": 1,         # Photometry for the 5-sigma detections
                "prv_candidates.magpsf": 1,
                "prv_candidates.sigmapsf": 1,
                "prv_candidates.band": 1,
                "prv_candidates.snr": 1,
                "prv_candidates.ra": 1,
                "prv_candidates.dec": 1,

                "fp_hists.band": 1,             # Forced photometry
                "fp_hists.jd": 1,
                "fp_hists.magpsf": 1,
                "fp_hists.sigmapsf": 1,
                "fp_hists.snr_psf": 1
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


# =========================
# SUMMARY HELPERS
# =========================


def make_run_dir(output_dir, start_time, end_time):
    savetime = Time.now().iso.replace(':', '-').replace(' ', '_')
    range_string = f"{start_time.replace(':','-')}_to_{end_time.replace(':','-')}"
    run_id = f"run_{savetime}_range_{range_string}"
    run_dir = os.path.join(output_dir, run_id)
    os.makedirs(run_dir, exist_ok=False)
    return run_dir



def as_list(value):
    if isinstance(value, list):
        return value
    return []



def finite_values(values):
    clean = []
    for value in values:
        if isinstance(value, (int, float)):
            clean.append(value)
    return clean



def first_non_null(points, key):
    for point in points:
        value = point.get(key)
        if value is not None:
            return value
    return None



def band_label(value):
    if value is None:
        return "unknown"

    label = str(value)
    safe = []
    for char in label:
        if char.isalnum():
            safe.append(char)
        else:
            safe.append("_")

    safe_label = "".join(safe).strip("_")
    return safe_label if safe_label else "unknown"



def min_or_none(values):
    values = finite_values(values)
    return min(values) if values else None



def max_or_none(values):
    values = finite_values(values)
    return max(values) if values else None



def build_candidate_summary(candidate_jsons):
    candidate_list = []

    for candidate in candidate_jsons:
        cand = candidate.get("candidate", {})

        candidate_list.append({
            "candid": candidate.get("_id"),
            "objectid": candidate.get("objectId"),
            "jd": cand.get("jd"),
            "magpsf": cand.get("magpsf"),
            "snr_psf": cand.get("snr_psf", cand.get("snr")),
        })

    return candidate_list



def write_candidate_summary(candidate_list, output_csv):
    fieldnames = ["candid", "objectid", "jd", "magpsf", "snr_psf"]

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(candidate_list)



def build_object_summary(object_jsons, candidate_jsons):
    n_alerts_by_object = Counter(candidate.get("objectId") for candidate in candidate_jsons)
    rows = []
    all_band_columns = set()

    for obj in object_jsons:
        object_id = obj.get("objectId", obj.get("_id"))

        prv_candidates = as_list(obj.get("prv_candidates"))
        fp_hists = as_list(obj.get("fp_hists"))
        all_points = prv_candidates + fp_hists

        prv_jds = finite_values(point.get("jd") for point in prv_candidates)
        fp_jds = finite_values(point.get("jd") for point in fp_hists)
        all_jds = prv_jds + fp_jds

        prv_mags = finite_values(point.get("magpsf") for point in prv_candidates)
        fp_mags = finite_values(point.get("magpsf") for point in fp_hists)
        all_mags = prv_mags + fp_mags

        prv_snrs = finite_values(point.get("snr") for point in prv_candidates)
        fp_snrs = finite_values(point.get("snr_psf") for point in fp_hists)
        all_snrs = prv_snrs + fp_snrs

        row = {
            "objectId": object_id,
            "n_alerts_in_run": n_alerts_by_object.get(object_id, 0),
            "n_prv_candidates": len(prv_candidates),
            "n_fp_hists": len(fp_hists),
            "n_phot_points_total": len(all_points),
            "t0_jd": min_or_none(all_jds),
            "t_last_jd": max_or_none(all_jds),
            "baseline_days": None,
            "t0_detection_jd": min_or_none(prv_jds),
            "t_last_detection_jd": max_or_none(prv_jds),
            "t0_forced_jd": min_or_none(fp_jds),
            "t_last_forced_jd": max_or_none(fp_jds),
            "ra": first_non_null(prv_candidates, "ra"),
            "dec": first_non_null(prv_candidates, "dec"),
            "brightest_magpsf": min_or_none(all_mags),
            "faintest_magpsf": max_or_none(all_mags),
            "max_snr": max_or_none(all_snrs),
        }

        if row["t0_jd"] is not None and row["t_last_jd"] is not None:
            row["baseline_days"] = row["t_last_jd"] - row["t0_jd"]

        prv_counts_by_band = Counter(band_label(point.get("band")) for point in prv_candidates)
        fp_counts_by_band = Counter(band_label(point.get("band")) for point in fp_hists)
        bands = sorted(set(prv_counts_by_band) | set(fp_counts_by_band))

        for band in bands:
            prv_col = f"n_prv_band_{band}"
            fp_col = f"n_fp_band_{band}"
            total_col = f"n_total_band_{band}"

            row[prv_col] = prv_counts_by_band.get(band, 0)
            row[fp_col] = fp_counts_by_band.get(band, 0)
            row[total_col] = row[prv_col] + row[fp_col]

            all_band_columns.update([prv_col, fp_col, total_col])

        rows.append(row)

    fixed_columns = [
        "objectId",
        "n_alerts_in_run",
        "n_prv_candidates",
        "n_fp_hists",
        "n_phot_points_total",
        "t0_jd",
        "t_last_jd",
        "baseline_days",
        "t0_detection_jd",
        "t_last_detection_jd",
        "t0_forced_jd",
        "t_last_forced_jd",
        "ra",
        "dec",
        "brightest_magpsf",
        "faintest_magpsf",
        "max_snr",
    ]

    fieldnames = fixed_columns + sorted(all_band_columns)
    return rows, fieldnames



def write_object_summary(rows, fieldnames, output_csv):
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# =========================
# MAIN
# =========================



def main(output_dir, filter_file, start_time, end_time):

    print("Auth...")
    token, _ = get_token()

    start_dt = Time(start_time, scale='utc').to_datetime()
    end_dt = Time(end_time, scale='utc').to_datetime()

    all_alerts = []

    print(f"Loading filter from: {filter_file}")
    pipeline = load_pipeline(filter_file)
    print(f"Loaded pipeline with {len(pipeline)} stages")

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
        if len(alerts) >= 10000:
            print("  WARNING: fetched 10000 alerts; the BOOM query limit may have truncated this batch.")

        all_alerts.extend(alerts)

        current = next_day

    print(f"\n                  Total alerts collected: {len(all_alerts)}")

    id_list = [a["_id"] for a in all_alerts]
    objectid_list = list(dict.fromkeys(a["objectId"] for a in all_alerts))

    print("===================================================================")
    print(f"                => Unique candidate IDs: {len(id_list)}")
    print(f"                => Unique object IDs: {len(objectid_list)}")
    print("===================================================================\n")

    run_dir = make_run_dir(output_dir, start_time, end_time)
    print(f"Saving this run to: {run_dir}\n")

    candidate_jsons = get_candidates(token, id_list)
    object_jsons = get_objects(token, objectid_list)
    print("")

    output_candidates_json = os.path.join(run_dir, "cands.json")
    output_objects_json = os.path.join(run_dir, "objs.json")
    output_candidates_summary = os.path.join(run_dir, "cands_sum.csv")
    output_objects_summary = os.path.join(run_dir, "objs_sum.csv")

    # =========================
    # BUILD AND SAVE SUMMARIES
    # =========================

    candidate_list = build_candidate_summary(candidate_jsons)
    write_candidate_summary(candidate_list, output_candidates_summary)
    print(f"Saved candidate summary to {output_candidates_summary}")

    object_summary_rows, object_summary_fieldnames = build_object_summary(object_jsons, candidate_jsons)
    write_object_summary(object_summary_rows, object_summary_fieldnames, output_objects_summary)
    print(f"Saved object summary to {output_objects_summary}")

    # =========================
    # SAVE JSON
    # =========================

    with open(output_candidates_json, "w") as f:
        json.dump(candidate_jsons, f, indent=2)
    print(f"Saved candidates to {output_candidates_json}")

    with open(output_objects_json, "w") as f:
        json.dump(object_jsons, f, indent=2)
    print(f"Saved objects to {output_objects_json}")

    return {
        "run_dir": run_dir,
        "candidates_json": output_candidates_json,
        "objects_json": output_objects_json,
        "candidates_summary": output_candidates_summary,
        "objects_summary": output_objects_summary,
        "n_alerts": len(all_alerts),
        "n_candidates": len(id_list),
        "n_objects": len(objectid_list),
    }



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run BOOM filter")
    parser.add_argument("--output_dir", default="boom_filter_outputs", help="Base output directory. A new run subfolder is created inside it.")
    parser.add_argument("--filter_file", required=True, help="Path to MongoDB pipeline JSON file")
    parser.add_argument("--start_time", required=True, help="ISO time")
    parser.add_argument("--end_time", required=True, help="ISO time")

    args = parser.parse_args()

    main(
        args.output_dir,
        args.filter_file,
        args.start_time,
        args.end_time
    )
