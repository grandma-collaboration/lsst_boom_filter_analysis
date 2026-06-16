#!/usr/bin/env python3
"""
Write a markdown summary for one BOOM pipeline run.
"""

from __future__ import annotations

import argparse
import csv
import os


def load_object_summary(objects_summary_csv):
    with open(objects_summary_csv, "r", newline="") as f:
        return list(csv.DictReader(f))


def load_delta_metrics(delta_metrics_csv):
    with open(delta_metrics_csv, "r", newline="") as f:
        return list(csv.DictReader(f))


def float_or_none(value):
    if value in (None, ""):
        return None
    return float(value)


def int_or_zero(value):
    if value in (None, ""):
        return 0
    return int(float(value))


def build_at_least_counts(object_rows):
    rows_with_counts = [
        {
            "objectId": row.get("objectId", ""),
            "n_alerts_in_run": int_or_zero(row.get("n_alerts_in_run")),
        }
        for row in object_rows
    ]
    max_count = max((row["n_alerts_in_run"] for row in rows_with_counts), default=0)
    n_total_objects = len(object_rows)

    at_least_rows = []
    for n_points in range(1, max_count + 1):
        matching_rows = [
            row for row in rows_with_counts
            if row["n_alerts_in_run"] >= n_points
        ]
        n_objects = len(matching_rows)
        proportion = 100 * n_objects / n_total_objects if n_total_objects > 0 else 0
        at_least_rows.append({
            "n_points": n_points,
            "n_objects": n_objects,
            "proportion_percent": proportion,
            "first_object_ids": [
                row["objectId"] for row in matching_rows[:5]
            ],
        })

    return at_least_rows


def first_relative_time_range_bin_objects(delta_metric_rows):
    matching_object_ids = []
    for row in delta_metric_rows:
        baseline_days = float_or_none(row.get("baseline_days"))
        relative_time_range = float_or_none(row.get("delta"))
        if baseline_days is None or relative_time_range is None:
            continue
        if baseline_days < 30 and 0 <= relative_time_range < 0.05:
            matching_object_ids.append(row.get("objectId", ""))

    return matching_object_ids


def markdown_table(headers, rows):
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def build_summary_markdown(n_candidates, n_objects, at_least_rows, first_bin_object_ids=None):
    sections = [
        "# Query Summary",
        "",
        "## Retrieved",
        "",
        markdown_table(
            ["Quantity", "Count"],
            [
                ["Candidates", n_candidates],
                ["Objects", n_objects],
            ],
        ),
        "",
        "## Objects With At Least N Filter-Passing Points",
        "",
    ]

    if at_least_rows:
        sections.append(markdown_table(
            [
                "At least N points passed filter",
                "Number of objects",
                "Proportion of total objects",
                "First 20 objectIds",
            ],
            [
                [
                    row["n_points"],
                    row["n_objects"],
                    f"{row['proportion_percent']:.2f}%",
                    ", ".join(row["first_object_ids"]),
                ]
                for row in at_least_rows
            ],
        ))
    else:
        sections.append("No objects with filter-passing points were found.")

    if first_bin_object_ids is not None:
        sections.extend([
            "",
            "## Objects In First Relative Time Range Bin",
            "",
            "Objects listed here have `baseline_days < 30` and `relative time range of interest < 0.05`.",
            "",
            markdown_table(
                ["Quantity", "Value"],
                [
                    ["Matching objects", len(first_bin_object_ids)],
                    ["First 40 objectIds", ", ".join(first_bin_object_ids[:40])],
                ],
            ),
        ])

    sections.append("")

    return "\n".join(sections)


def main(
    objects_summary_csv,
    output_dir,
    n_candidates,
    n_objects,
    output_filename="query_summary.md",
    delta_metrics_csv=None,
):
    os.makedirs(output_dir, exist_ok=True)

    object_rows = load_object_summary(objects_summary_csv)
    at_least_rows = build_at_least_counts(object_rows)
    first_bin_object_ids = None
    if delta_metrics_csv is not None:
        first_bin_object_ids = first_relative_time_range_bin_objects(
            load_delta_metrics(delta_metrics_csv)
        )

    markdown = build_summary_markdown(
        n_candidates=n_candidates,
        n_objects=n_objects,
        at_least_rows=at_least_rows,
        first_bin_object_ids=first_bin_object_ids,
    )

    output_path = os.path.join(output_dir, output_filename)
    with open(output_path, "w") as f:
        f.write(markdown)

    print(f"Saved query summary to {output_path}")
    return {
        "query_summary_markdown": output_path,
        "n_objects_at_least_one_passing_point": at_least_rows[0]["n_objects"] if at_least_rows else 0,
        "max_filter_passing_points_per_object": at_least_rows[-1]["n_points"] if at_least_rows else 0,
        "n_objects_first_relative_time_range_bin": len(first_bin_object_ids) if first_bin_object_ids is not None else None,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write a markdown summary for one BOOM pipeline run.")
    parser.add_argument("--objects_summary_csv", required=True, help="Path to objs_sum.csv")
    parser.add_argument("--output_dir", required=True, help="Run directory where query_summary.md is written")
    parser.add_argument("--n_candidates", type=int, required=True, help="Total candidates retrieved")
    parser.add_argument("--n_objects", type=int, required=True, help="Total objects retrieved")
    parser.add_argument("--output_filename", default="query_summary.md", help="Summary markdown filename")
    parser.add_argument("--delta_metrics_csv", default=None, help="Optional path to plots/delta_metrics.csv")
    args = parser.parse_args()

    main(
        args.objects_summary_csv,
        args.output_dir,
        n_candidates=args.n_candidates,
        n_objects=args.n_objects,
        output_filename=args.output_filename,
        delta_metrics_csv=args.delta_metrics_csv,
    )
