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
    band_rates.json
"""

# =========================
# CONFIG
# =========================

BOOM_URL = "https://api.kaboom.caltech.edu"

# Most filter-passing alerts BOOM returns for one 1-day filter query
FILTER_QUERY_LIMIT = 10000

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
        # "object_ids": ["170635513149522077"],
        "start_jd": start_jd,
        "end_jd": end_jd,
        "survey": "LSST",
        "permissions": {"LSST": []},
        "limit": FILTER_QUERY_LIMIT,
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

                "prv_candidates.diaSourceId": 1,         # Photometry for the 5-sigma detections
                "prv_candidates.jd": 1,
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


def get_last_point_ids_by_band(object_jsons):
    """
    For each object, find the diaSourceId of its most recent (max jd) LSST
    detection in every band. Each LSST detection is itself an alert document
    in LSST_alerts, and its diaSourceId doubles as that document's _id, so
    this identifies exactly which alert to query for the rate as of that
    specific point -- not the rate of the last point that happened to pass a
    filter, which may be an earlier, now-stale fit.
    """
    last_point_by_object_band = {}
    for obj in object_jsons:
        obj_id = obj.get("objectId", obj.get("_id"))
        best = {}  # band -> (jd, diaSourceId)
        for point in as_list(obj.get("prv_candidates")):
            band = point.get("band")
            jd = point.get("jd")
            dia_source_id = point.get("diaSourceId")
            if band is None or jd is None or dia_source_id is None:
                continue
            current = best.get(band)
            if current is None or jd > current[0]:
                best[band] = (jd, dia_source_id)
        last_point_by_object_band[obj_id] = {
            band: dia_source_id for band, (jd, dia_source_id) in best.items()
        }
    return last_point_by_object_band


def get_band_rates(token, object_jsons, batch_size=1000):
    """
    Fetch BOOM's rising/fading rate for the most recent LSST detection in
    each band of each object, regardless of whether that detection ever
    passed a filter.

    LSST_alerts_aux (the source of prv_candidates/objs.json) carries only bare
    photometry: no properties, no photstats, on the document or on any
    individual point. The rate only exists on LSST_alerts documents -- one per
    detection, each with a properties.photstats snapshot as of that
    detection -- so this queries LSST_alerts directly by diaSourceId for
    exactly the points identified by get_last_point_ids_by_band.
    """
    last_point_by_object_band = get_last_point_ids_by_band(object_jsons)

    dia_source_id_to_keys = {}  # diaSourceId -> list of (obj_id, band)
    for obj_id, bands in last_point_by_object_band.items():
        for band, dia_source_id in bands.items():
            dia_source_id_to_keys.setdefault(dia_source_id, []).append((obj_id, band))

    dia_source_ids = list(dia_source_id_to_keys.keys())

    url = f"{BOOM_URL}/queries/find"
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }

    rates_by_object = {}
    for i in range(0, len(dia_source_ids), batch_size):
        batch_ids = dia_source_ids[i:i + batch_size]

        payload = {
            "catalog_name": "LSST_alerts",
            "filter": {"_id": {"$in": batch_ids}},
            "projection": {
                "objectId": 1,
                "candidate.band": 1,
                "properties.photstats": 1,
            },
        }

        print(f"Querying band-rate batch {i // batch_size + 1} ({len(batch_ids)} ids)")

        response = requests.post(url, json=payload, headers=headers)

        if response.status_code != 200:
            print(response.status_code)
            print(response.text)
            response.raise_for_status()

        for alert in response.json()["data"]:
            dia_source_id = alert.get("_id")
            photstats = (alert.get("properties") or {}).get("photstats") or {}
            for obj_id, band in dia_source_id_to_keys.get(dia_source_id, []):
                band_stats = photstats.get(band) or {}
                rising = band_stats.get("rising") or {}
                fading = band_stats.get("fading") or {}
                rates_by_object.setdefault(obj_id, {})[band] = {
                    "rising_rate": rising.get("rate"),
                    "fading_rate": fading.get("rate"),
                }

    return rates_by_object


# =========================
# TOTAL COUNTS
# =========================


def count_documents(token, catalog_name, filter_dict):
    """Number of documents in catalog_name matching filter_dict (BOOM /queries/count)."""
    url = f"{BOOM_URL}/queries/count"

    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }

    payload = {
        "catalog_name": catalog_name,
        "filter": filter_dict,
    }

    response = requests.post(url, json=payload, headers=headers)

    if response.status_code != 200:
        print(response.status_code)
        print(response.text)
        response.raise_for_status()

    return response.json()["data"]


def jd_range(start_time, end_time):
    """Mongo range on a Julian Date for the half-open UTC interval [start_time, end_time)."""
    return {
        "$gte": Time(start_time, scale='utc').jd,
        "$lt": Time(end_time, scale='utc').jd,
    }


def build_total_count_filters(start_time, end_time, exclude_start_time=None, exclude_end_time=None):
    """
    Filters selecting every LSST alert and every LSST object in the queried
    range -- the populations the filter-passing alerts and objects are drawn
    from, so they can be turned into pass fractions:
      - alerts (LSST_alerts): detection time in [start, end), minus the
        alerts inside the excluded range;
      - objects (LSST_alerts_aux): at least one detection in [start, end),
        and no point at all -- detection or forced photometry -- in the
        excluded range, the same rule drop_objects_with_excluded_points
        applies to the filter-passing objects.
    $elemMatch keeps each range on a single point, so an object with points
    on both sides of the excluded range, but none inside it, is kept.
    """
    in_range = jd_range(start_time, end_time)
    alert_conditions = [{"candidate.jd": in_range}]
    object_conditions = [{"prv_candidates": {"$elemMatch": {"jd": in_range}}}]

    if exclude_start_time and exclude_end_time:
        excluded = jd_range(exclude_start_time, exclude_end_time)
        alert_conditions.append({"$nor": [{"candidate.jd": excluded}]})
        object_conditions.append({"$nor": [
            {"prv_candidates": {"$elemMatch": {"jd": excluded}}},
            {"fp_hists": {"$elemMatch": {"jd": excluded}}},
        ]})

    return {"$and": alert_conditions}, {"$and": object_conditions}


def pass_fraction(n_passing, n_total):
    return n_passing / n_total if n_total else None


def format_pass_fraction(n_passing, n_total):
    fraction = pass_fraction(n_passing, n_total)
    percent = f"{100 * fraction:.3g}%" if fraction is not None else "n/a"
    return f"{n_passing:,} / {n_total:,} ({percent})"


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


def build_query_intervals(start_dt, end_dt, exclude_start_dt=None, exclude_end_dt=None):
    if exclude_start_dt is None and exclude_end_dt is None:
        return [(start_dt, end_dt)]

    if exclude_start_dt is None or exclude_end_dt is None:
        raise ValueError("exclude_start_time and exclude_end_time must be provided together")

    if exclude_end_dt <= exclude_start_dt:
        raise ValueError("exclude_end_time must be after exclude_start_time")

    if exclude_end_dt <= start_dt or exclude_start_dt >= end_dt:
        return [(start_dt, end_dt)]

    intervals = []
    if exclude_start_dt > start_dt:
        intervals.append((start_dt, min(exclude_start_dt, end_dt)))
    if exclude_end_dt < end_dt:
        intervals.append((max(exclude_end_dt, start_dt), end_dt))

    if not intervals:
        raise ValueError("Excluded time range removes the entire query interval")

    return intervals


def interval_to_dict(start_dt, end_dt):
    return {
        "start_time": start_dt.isoformat(),
        "end_time": end_dt.isoformat(),
    }



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



def make_exclude_jd(exclude_start_time, exclude_end_time):
    if not (exclude_start_time and exclude_end_time):
        return None
    return (
        Time(exclude_start_time, scale='utc').jd,
        Time(exclude_end_time, scale='utc').jd,
    )


def object_has_excluded_point(obj, exclude_jd):
    """
    True if this object has any point, in any band, in either prv_candidates
    (LSST detections) or fp_hists (forced photometry), inside the half-open
    excluded interval [start_jd, end_jd).
    """
    if exclude_jd is None:
        return False
    start_jd, end_jd = exclude_jd
    for key in ("prv_candidates", "fp_hists"):
        for point in as_list(obj.get(key)):
            jd = point.get("jd")
            if isinstance(jd, (int, float)) and start_jd <= jd < end_jd:
                return True
    return False


def drop_objects_with_excluded_points(object_jsons, exclude_jd):
    """
    Drop every object that has any point in the excluded interval, instead of
    trimming individual points. An object is either kept with its full,
    untouched history, or left out of the run entirely -- so there is no
    partial history left behind, and no baseline to reconcile downstream.
    """
    if exclude_jd is None:
        return object_jsons, set()

    kept = []
    dropped_object_ids = set()
    for obj in object_jsons:
        object_id = obj.get("objectId", obj.get("_id"))
        if object_has_excluded_point(obj, exclude_jd):
            dropped_object_ids.add(object_id)
        else:
            kept.append(obj)
    return kept, dropped_object_ids



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



def main(output_dir, filter_file, start_time, end_time, exclude_start_time=None, exclude_end_time=None, compute_pass_rates=False):

    print("Auth...")
    token, _ = get_token()

    start_dt = Time(start_time, scale='utc').to_datetime()
    end_dt = Time(end_time, scale='utc').to_datetime()
    exclude_start_dt = Time(exclude_start_time, scale='utc').to_datetime() if exclude_start_time else None
    exclude_end_dt = Time(exclude_end_time, scale='utc').to_datetime() if exclude_end_time else None
    query_intervals = build_query_intervals(
        start_dt,
        end_dt,
        exclude_start_dt=exclude_start_dt,
        exclude_end_dt=exclude_end_dt,
    )

    all_alerts = []

    print(f"Loading filter from: {filter_file}")
    pipeline = load_pipeline(filter_file)
    print(f"Loaded pipeline with {len(pipeline)} stages")

    print("Running LSST filter in 1-day batches...\n")
    if exclude_start_time and exclude_end_time:
        print(f"Excluding time range: {exclude_start_time} to {exclude_end_time}")

    for interval_start, interval_end in query_intervals:
        print(f"Included query interval: {interval_start.isoformat()} to {interval_end.isoformat()}")
        current = interval_start

        while current < interval_end:
            next_day = min(current + timedelta(days=1), interval_end)

            start_iso = current.isoformat()
            end_iso = next_day.isoformat()

            print(f"Querying: {start_iso} to {end_iso}")

            result = run_filter(token, pipeline, start_iso, end_iso)
            alerts = result.get("data", {}).get("results", [])

            print(f"  fetched {len(alerts)} alerts")
            if len(alerts) >= FILTER_QUERY_LIMIT:
                print(f"  WARNING: fetched {FILTER_QUERY_LIMIT} alerts; the BOOM query limit may have truncated this batch.")

            all_alerts.extend(alerts)

            current = next_day

    print(f"\n                  Total alerts collected: {len(all_alerts)}")

    id_list = [a["_id"] for a in all_alerts]
    objectid_list = list(dict.fromkeys(a["objectId"] for a in all_alerts))

    run_dir = make_run_dir(output_dir, start_time, end_time)
    print(f"Saving this run to: {run_dir}\n")

    candidate_jsons = get_candidates(token, id_list)
    object_jsons = get_objects(token, objectid_list)
    print("")

    # An object with any point (any band, detection or forced photometry) in
    # the excluded range is dropped entirely, along with its candidates, so
    # this run's outputs behave as if such objects were never retrieved.
    exclude_jd = make_exclude_jd(exclude_start_time, exclude_end_time)
    object_jsons, dropped_object_ids = drop_objects_with_excluded_points(object_jsons, exclude_jd)
    if dropped_object_ids:
        candidate_jsons = [c for c in candidate_jsons if c.get("objectId") not in dropped_object_ids]
        print(
            f"Dropped {len(dropped_object_ids)} object(s) with photometry in the excluded range "
            f"({exclude_start_time} to {exclude_end_time}); "
            f"{len(object_jsons)} object(s) and {len(candidate_jsons)} candidate(s) remain\n"
        )

    # =========================
    # PASS FRACTIONS
    # =========================

    # Alerts: the excluded range only removes the alerts inside it, which the
    # filter never queried, so every filter-passing alert counts.
    # Objects: objects with points in the excluded range are not in the total,
    # so only the ones left after the drop above count.
    n_alerts_passing = len(dict.fromkeys(id_list))
    n_objects_passing = len(object_jsons)

    filter_pass_rates = None
    if compute_pass_rates:
        # Every LSST alert and object in the queried range, with the excluded
        # range applied the same way as to the filter-passing ones. The alert
        # count is indexed; the object count scans all of LSST_alerts_aux (no
        # index on prv_candidates.jd or fp_hists.jd), which takes about a minute.
        alerts_total_filter, objects_total_filter = build_total_count_filters(
            start_time, end_time, exclude_start_time, exclude_end_time
        )
        print("Counting every LSST alert in the queried range...")
        n_alerts_total = count_documents(token, "LSST_alerts", alerts_total_filter)
        print("Counting every LSST object in the queried range...")
        n_objects_total = count_documents(token, "LSST_alerts_aux", objects_total_filter)

        filter_pass_rates = {
            "n_alerts_passing": n_alerts_passing,
            "n_alerts_total": n_alerts_total,
            "alert_pass_fraction": pass_fraction(n_alerts_passing, n_alerts_total),
            "n_objects_passing": n_objects_passing,
            "n_objects_total": n_objects_total,
            "object_pass_fraction": pass_fraction(n_objects_passing, n_objects_total),
            "alerts_total_count_query": {"catalog_name": "LSST_alerts", "filter": alerts_total_filter},
            "objects_total_count_query": {"catalog_name": "LSST_alerts_aux", "filter": objects_total_filter},
        }

    print("\n===================================================================")
    if filter_pass_rates is not None:
        print(f"                => Filter-passing alerts:  {format_pass_fraction(n_alerts_passing, n_alerts_total)}")
        print(f"                => Filter-passing objects: {format_pass_fraction(n_objects_passing, n_objects_total)}")
    else:
        print(f"                => Filter-passing alerts:  {n_alerts_passing:,}")
        print(f"                => Filter-passing objects: {n_objects_passing:,}")
    print("===================================================================\n")

    print("Querying the latest rising/fading rate per band for each object...")
    band_rates = get_band_rates(token, object_jsons)

    output_candidates_json = os.path.join(run_dir, "cands.json")
    output_objects_json = os.path.join(run_dir, "objs.json")
    output_candidates_summary = os.path.join(run_dir, "cands_sum.csv")
    output_objects_summary = os.path.join(run_dir, "objs_sum.csv")
    output_band_rates_json = os.path.join(run_dir, "band_rates.json")

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

    with open(output_band_rates_json, "w") as f:
        json.dump(band_rates, f, indent=2)
    print(f"Saved band rates to {output_band_rates_json}")

    return {
        "run_dir": run_dir,
        "candidates_json": output_candidates_json,
        "objects_json": output_objects_json,
        "candidates_summary": output_candidates_summary,
        "objects_summary": output_objects_summary,
        "band_rates_json": output_band_rates_json,
        "n_alerts": len(all_alerts),
        "n_candidates": len(candidate_jsons),
        "n_objects": len(object_jsons),
        "n_objects_dropped_excluded": len(dropped_object_ids),
        "filter_pass_rates": filter_pass_rates,
        "exclude_start_time": exclude_start_time,
        "exclude_end_time": exclude_end_time,
        "query_intervals": [
            interval_to_dict(interval_start, interval_end)
            for interval_start, interval_end in query_intervals
        ],
    }



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run BOOM filter")
    parser.add_argument("--output_dir", default="boom_filter_outputs", help="Base output directory. A new run subfolder is created inside it.")
    parser.add_argument("--filter_file", required=True, help="Path to MongoDB pipeline JSON file")
    parser.add_argument("--start_time", required=True, help="ISO time")
    parser.add_argument("--end_time", required=True, help="ISO time")
    parser.add_argument("--exclude_start_time", default=None, help="Optional ISO start time of a range to exclude: it is not queried, and objects with any point in it are dropped")
    parser.add_argument("--exclude_end_time", default=None, help="Optional ISO end time of the excluded range")
    parser.add_argument("--compute_pass_rates", action="store_true", help="Also count every LSST alert/object in the range and report the filter pass fractions (adds about 2 minutes)")

    args = parser.parse_args()

    main(
        args.output_dir,
        args.filter_file,
        args.start_time,
        args.end_time,
        exclude_start_time=args.exclude_start_time,
        exclude_end_time=args.exclude_end_time,
        compute_pass_rates=args.compute_pass_rates,
    )
