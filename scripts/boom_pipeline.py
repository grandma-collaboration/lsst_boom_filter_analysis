#!/usr/bin/env python3
"""
Run the complete BOOM light-curve pipeline in one command.

The pipeline:
  1. Runs a BOOM filter over the requested time range.
  2. Saves candidates/objects JSON and CSV summaries in one run directory.
  3. Plots LSST detections, LSST forced photometry, and optional ZTF crossmatch
     photometry into the same run directory.

Example:
  python3 boom_pipeline.py \
    --output_dir ../boom_filter_outputs \
    --filter_file ../filters/my_filter.json \
    --start_time 2026-03-01T12:00:00 \
    --end_time 2026-05-13T12:00:00
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from astropy.time import Time


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


def run_pipeline(args: argparse.Namespace) -> dict[str, Any]:
    if not os.environ.get("BOOM_PASSWORD"):
        raise RuntimeError("BOOM_PASSWORD environment variable not set")
    
    if not os.environ.get("BOOM_USERNAME"):
        raise RuntimeError("BOOM_USERNAME environment variable not set")

    filter_file = resolve_existing_file(args.filter_file, "Filter file")
    validate_time_range(args.start_time, args.end_time)
    if args.max_objects is not None and args.max_objects < 1:
        raise ValueError("max_objects must be a positive integer")

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Imported after validation so command-line errors stay clear and fast.
    import boom_filter_per_run
    import plot_lc

    print("\n========== BOOM filter ==========")
    filter_result = boom_filter_per_run.main(
        str(output_dir),
        str(filter_file),
        args.start_time,
        args.end_time,
    )

    run_dir = Path(filter_result["run_dir"])
    plots_subdir = validate_relative_subdir(args.plots_subdir, "plots_subdir")

    print("\n========== Light-curve plots ==========")
    plot_result = plot_lc.main(
        filter_result["objects_json"],
        filter_result["candidates_json"],
        str(run_dir),
        include_ztf=not args.no_ztf,
        include_fp=not args.no_forced_photometry,
        output_subdir=plots_subdir,
        max_objects=args.max_objects,
    )

    manifest = {
        "created_at_utc": iso_utc_now(),
        "status": "complete",
        "command": {
            "filter_file": str(filter_file),
            "start_time": args.start_time,
            "end_time": args.end_time,
            "include_ztf": not args.no_ztf,
            "include_forced_photometry": not args.no_forced_photometry,
            "max_objects": args.max_objects,
        },
        "filter_result": filter_result,
        "plot_result": plot_result,
        "artifacts": {
            "run_dir": str(run_dir),
            "candidates_json": filter_result["candidates_json"],
            "objects_json": filter_result["objects_json"],
            "candidates_summary_csv": filter_result["candidates_summary"],
            "objects_summary_csv": filter_result["objects_summary"],
            "plots_dir": plot_result["plots_dir"],
            "manifest": str(run_dir / "pipeline_manifest.json"),
        },
    }
    manifest_path = write_manifest(run_dir, manifest)
    manifest["artifacts"]["manifest"] = str(manifest_path)

    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run BOOM filtering and complete light-curve plotting in one command."
    )
    parser.add_argument(
        "--output_dir",
        default="boom_filter_outputs",
        help="Base output directory. A timestamped run folder is created inside it.",
    )
    parser.add_argument("--filter_file", required=True, help="Path to MongoDB pipeline JSON file")
    parser.add_argument("--start_time", required=True, help="UTC ISO start time")
    parser.add_argument("--end_time", required=True, help="UTC ISO end time")
    parser.add_argument(
        "--plots_subdir",
        default="plots",
        help="Subdirectory inside the run folder where plots are written.",
    )
    parser.add_argument(
        "--no_ztf",
        action="store_true",
        help="Disable ZTF cone-search crossmatch photometry in the plots.",
    )
    parser.add_argument(
        "--no_forced_photometry",
        action="store_true",
        help="Disable LSST forced photometry in the plots.",
    )
    parser.add_argument(
        "--max_objects",
        type=int,
        default=None,
        help="Optional cap on the number of objects to plot, useful for smoke tests.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        manifest = run_pipeline(args)
    except Exception as exc:
        print(f"\nPipeline failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    artifacts = manifest["artifacts"]
    print("\n========== Pipeline complete ==========")
    print(f"Run directory:          {artifacts['run_dir']}")
    print(f"Candidate summary CSV:  {artifacts['candidates_summary_csv']}")
    print(f"Candidate JSON:         {artifacts['candidates_json']}")
    print(f"Object JSON:            {artifacts['objects_json']}")
    print(f"Object summary CSV:     {artifacts['objects_summary_csv']}")
    print(f"Plots directory:        {artifacts['plots_dir']}")
    print(f"Manifest:               {artifacts['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
