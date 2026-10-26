#!/usr/bin/env python3
"""
Run the complete BOOM light-curve pipeline in one command.

The pipeline:
  1. Runs a BOOM filter over the requested time range.
  2. Saves candidates/objects JSON and CSV summaries in one run directory.
  3. Plots LSST detections, LSST forced photometry, and optional ZTF crossmatch
     photometry into the same run directory.

There is no command-line interface for this entry point. Every parameter
lives in inputs.py, which is the only file to edit before a run:

  1. Edit scripts/inputs.py.
  2. Run:
       python3 scripts/boom_pipeline.py

(The individual step scripts this pipeline calls -- boom_filter_per_run.py,
plot_lc.py, population_plots.py, run_summary.py -- can still be run on their
own from the command line with their own arguments; see their --help.)
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from astropy.time import Time

import inputs


def iso_utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def resolve_existing_file(path: str, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist or is not a file: {resolved}")
    return resolved


def validate_time_range(start_time: str, end_time: str) -> None:
    start = Time(start_time, scale="utc")
    end = Time(end_time, scale="utc")
    if end <= start:
        raise ValueError(f"end_time must be after start_time ({start_time} >= {end_time})")


def validate_exclude_time_range(
    start_time: str,
    end_time: str,
    exclude_start_time: str | None,
    exclude_end_time: str | None,
) -> None:
    if exclude_start_time is None and exclude_end_time is None:
        return

    if exclude_start_time is None or exclude_end_time is None:
        raise ValueError("EXCLUDE_START_TIME and EXCLUDE_END_TIME must be provided together")

    start = Time(start_time, scale="utc")
    end = Time(end_time, scale="utc")
    exclude_start = Time(exclude_start_time, scale="utc")
    exclude_end = Time(exclude_end_time, scale="utc")

    if exclude_end <= exclude_start:
        raise ValueError("EXCLUDE_END_TIME must be after EXCLUDE_START_TIME")

    if exclude_start <= start and exclude_end >= end:
        raise ValueError("Excluded time range removes the entire query interval")


def validate_relative_subdir(value: str, label: str) -> str:
    cleaned = value.strip().strip("/\\") or "plots"
    path = Path(cleaned)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"{label} must stay inside the run directory: {value}")
    return cleaned


def write_manifest(run_dir: str | Path, manifest: dict[str, Any]) -> Path:
    manifest_path = Path(run_dir) / "pipeline_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest_path


def run_pipeline() -> dict[str, Any]:
    if not os.environ.get("BOOM_PASSWORD"):
        raise RuntimeError("BOOM_PASSWORD environment variable not set")

    if not os.environ.get("BOOM_USERNAME"):
        raise RuntimeError("BOOM_USERNAME environment variable not set")

    filter_file = resolve_existing_file(inputs.FILTER_FILE, "FILTER_FILE")
    validate_time_range(inputs.START_TIME, inputs.END_TIME)
    validate_exclude_time_range(
        inputs.START_TIME,
        inputs.END_TIME,
        inputs.EXCLUDE_START_TIME,
        inputs.EXCLUDE_END_TIME,
    )
    if inputs.MAX_OBJECTS is not None and inputs.MAX_OBJECTS < 1:
        raise ValueError("MAX_OBJECTS must be a positive integer")

    output_dir = Path(inputs.OUTPUT_DIR).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Imported after validation so configuration errors stay clear and fast.
    import boom_filter_per_run
    import plot_lc
    import population_plots
    import run_summary

    print("\n========== BOOM filter ==========")
    filter_result = boom_filter_per_run.main(
        str(output_dir),
        str(filter_file),
        inputs.START_TIME,
        inputs.END_TIME,
        exclude_start_time=inputs.EXCLUDE_START_TIME,
        exclude_end_time=inputs.EXCLUDE_END_TIME,
        compute_pass_rates=inputs.COMPUTE_FILTER_PASS_RATES,
    )

    run_dir = Path(filter_result["run_dir"])
    plots_subdir = validate_relative_subdir(inputs.PLOTS_SUBDIR, "PLOTS_SUBDIR")

    if inputs.SKIP_INDIVIDUAL_PLOTS:
        print("\n========== LSST-only light-curve metrics ==========")
    else:
        print("\n========== Light-curve plots ==========")
    plot_result = plot_lc.main(
        filter_result["objects_json"],
        filter_result["candidates_json"],
        str(run_dir),
        include_ztf=inputs.INCLUDE_ZTF,
        include_fp=inputs.INCLUDE_FORCED_PHOTOMETRY,
        output_subdir=plots_subdir,
        max_objects=inputs.MAX_OBJECTS,
        make_plots=not inputs.SKIP_INDIVIDUAL_PLOTS,
        exclude_start_time=inputs.EXCLUDE_START_TIME,
        exclude_end_time=inputs.EXCLUDE_END_TIME,
        band_rates_json=filter_result["band_rates_json"],
    )

    print("\n========== Population plots ==========")
    population_plot_result = population_plots.main(
        plot_result["delta_metrics_csv"],
        str(run_dir),
        output_subdir="pop_plots",
    )

    print("\n========== Query summary ==========")
    summary_result = run_summary.main(
        filter_result["objects_summary"],
        str(run_dir),
        n_candidates=filter_result["n_candidates"],
        n_objects=filter_result["n_objects"],
        delta_metrics_csv=plot_result["delta_metrics_csv"],
    )

    manifest = {
        "created_at_utc": iso_utc_now(),
        "status": "complete",
        "command": {
            "filter_file": str(filter_file),
            "start_time": inputs.START_TIME,
            "end_time": inputs.END_TIME,
            "exclude_start_time": inputs.EXCLUDE_START_TIME,
            "exclude_end_time": inputs.EXCLUDE_END_TIME,
            "include_ztf": inputs.INCLUDE_ZTF,
            "include_forced_photometry": inputs.INCLUDE_FORCED_PHOTOMETRY,
            "max_objects": inputs.MAX_OBJECTS,
            "skip_individual_plots": inputs.SKIP_INDIVIDUAL_PLOTS,
            "compute_filter_pass_rates": inputs.COMPUTE_FILTER_PASS_RATES,
        },
        # Filter-passing alerts/objects vs. every LSST alert/object in the
        # queried range; None unless COMPUTE_FILTER_PASS_RATES is True
        # "filter_pass_rates": filter_result["filter_pass_rates"],
        "filter_result": filter_result,
        "plot_result": plot_result,
        "population_plot_result": population_plot_result,
        "summary_result": summary_result,
        "artifacts": {
            "run_dir": str(run_dir),
            "candidates_json": filter_result["candidates_json"],
            "objects_json": filter_result["objects_json"],
            "band_rates_json": filter_result["band_rates_json"],
            "candidates_summary_csv": filter_result["candidates_summary"],
            "objects_summary_csv": filter_result["objects_summary"],
            "plots_dir": plot_result["plots_dir"],
            "delta_metrics_csv": plot_result["delta_metrics_csv"],
            "pop_plots_dir": population_plot_result["pop_plots_dir"],
            "delta_distribution_plot": population_plot_result["delta_distribution_plot"],
            "highlighted_fraction_by_band_plot": population_plot_result["highlighted_fraction_by_band_plot"],
            "highlighted_fraction_global_plot": population_plot_result["highlighted_fraction_global_plot"],
            "lsst_points_by_band_plot": population_plot_result["lsst_points_by_band_plot"],
            "lsst_points_global_plot": population_plot_result["lsst_points_global_plot"],
            "query_summary_markdown": summary_result["query_summary_markdown"],
            "manifest": str(run_dir / "pipeline_manifest.json"),
        },
    }
    manifest_path = write_manifest(run_dir, manifest)
    manifest["artifacts"]["manifest"] = str(manifest_path)

    return manifest


def main() -> int:
    try:
        manifest = run_pipeline()
    except Exception as exc:
        print(f"\nPipeline failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    artifacts = manifest["artifacts"]
    print("\n========== Pipeline complete ==========")
    print(f"Run directory:          {artifacts['run_dir']}")
    print(f"Candidate summary CSV:  {artifacts['candidates_summary_csv']}")
    print(f"Candidate JSON:         {artifacts['candidates_json']}")
    print(f"Object JSON:            {artifacts['objects_json']}")
    print(f"Band rates JSON:        {artifacts['band_rates_json']}")
    print(f"Object summary CSV:     {artifacts['objects_summary_csv']}")
    print(f"Plots directory:        {artifacts['plots_dir']}")
    print(f"Delta metrics CSV:      {artifacts['delta_metrics_csv']}")
    print(f"Population plots dir:   {artifacts['pop_plots_dir']}")
    print(f"Time range distribution:{artifacts['delta_distribution_plot']}")
    print(f"Band fraction plot:     {artifacts['highlighted_fraction_by_band_plot']}")
    print(f"Global fraction plot:   {artifacts['highlighted_fraction_global_plot']}")
    print(f"Band LSST count plot:   {artifacts['lsst_points_by_band_plot']}")
    print(f"Global LSST count plot: {artifacts['lsst_points_global_plot']}")
    print(f"Query summary:          {artifacts['query_summary_markdown']}")
    print(f"Manifest:               {artifacts['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
