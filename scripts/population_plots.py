#!/usr/bin/env python3
"""
Create population-level plots for a BOOM pipeline run.

Currently this produces population distributions from the per-object metrics
written by plot_lc.py.
"""

from __future__ import annotations

import argparse
import csv
import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BANDS = ["u", "g", "r", "i", "z", "y"]
DELTA_MAX_BASELINE_DAYS = 30


def save_figure(fig, output_path):
    fig.savefig(output_path)
    pickle_path = os.path.splitext(output_path)[0] + ".pickle"
    with open(pickle_path, "wb") as f:
        pickle.dump(fig, f)
    return pickle_path


def load_delta_metrics(delta_metrics_csv):
    with open(delta_metrics_csv, "r", newline="") as f:
        return list(csv.DictReader(f))


def float_or_none(value):
    if value in (None, ""):
        return None
    return float(value)


def plot_delta_distribution(delta_metrics, output_dir):

    n_excluded_long_baseline = 0
    for row in delta_metrics:
        baseline_days = float_or_none(row.get("baseline_days"))
        if baseline_days is not None and baseline_days > DELTA_MAX_BASELINE_DAYS:
            n_excluded_long_baseline += 1

    n_singlepoint = 0
    n_both = 0

    for row in delta_metrics:
        if row.get("delta_status") == "single_point":
            n_singlepoint += 1
            if float_or_none(row.get("baseline_days")) > DELTA_MAX_BASELINE_DAYS:
                n_both += 1


    used_deltas = [
        float_or_none(row.get("delta"))
        for row in delta_metrics
        if row.get("delta_status") == "used"
        and float_or_none(row.get("baseline_days")) < DELTA_MAX_BASELINE_DAYS
    ]

    used_deltas = [value for value in used_deltas if value is not None]

    n_used = len(used_deltas)
    denominator = len(delta_metrics)
    used_fraction = n_used / denominator if denominator > 0 else 0

    fig, ax = plt.subplots(figsize=(8, 5))

    if used_deltas:
        bins = np.linspace(0, 1, 21)
        ax.hist(used_deltas, bins=bins, color="steelblue", edgecolor="white")
    else:
        ax.text(
            0.5,
            0.5,
            "No usable multi-point objects",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
        )

    if denominator > 0:
        proportion_text = f"Used objects: {n_used}/{denominator} = {used_fraction:.1%}"
    else:
        proportion_text = "Used objects: 0/0"
    ax.text(
        0.98,
        0.95,
        proportion_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )
    ax.text(
        0.98,
        0.85,
        f"Excluded baseline > {DELTA_MAX_BASELINE_DAYS} d: {n_excluded_long_baseline}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )
    ax.text(
        0.98,
        0.75,
        f"Excluded single-passing-point objects: {n_singlepoint}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )
    ax.text(
        0.98,
        0.65,
        f"Objects excluded for both conditions: {n_both}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )


    ax.set_title("Distribution of Relative Time Range of Interest")
    ax.set_xlabel("Relative time range of interest")
    ax.set_ylabel("Number of objects")
    ax.set_xlim(0, 1)
    ax.grid(axis="y", alpha=0.25)

    output_path = os.path.join(output_dir, "delta_distribution.png")
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close(fig)

    return {
        "delta_distribution_plot": output_path,
        "n_used_deltas": n_used,
        "n_single_point_deltas": n_singlepoint,
        "used_delta_fraction": used_fraction,
        "n_delta_excluded_long_baseline": n_excluded_long_baseline,
    }


def values_from_column(rows, column):
    values = []
    for row in rows:
        value = float_or_none(row.get(column))
        if value is not None:
            values.append(value)
    return values


def integer_bins(values):
    if not values:
        return np.arange(-0.5, 1.5, 1)
    upper = int(max(values))
    return np.arange(-0.5, upper + 1.5, 1)


def logarithmic_bins(values, n_bins=20):
    positive_values = [value for value in values if value > 0]
    if not positive_values:
        return None

    lower = min(positive_values)
    upper = max(positive_values)
    if lower == upper:
        lower = lower / np.sqrt(10)
        upper = upper * np.sqrt(10)

    return np.logspace(np.log10(lower), np.log10(upper), n_bins + 1)


def annotate_used_objects(ax, n_used, n_total):
    if n_total > 0:
        text = f"Used objects: {n_used}/{n_total} = {n_used / n_total:.1%}"
    else:
        text = "Used objects: 0/0"

    ax.text(
        0.98,
        0.95,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )


def plot_fraction_histogram(ax, values, title, n_total):
    if values:
        bins = np.linspace(0, 1, 21)
        ax.hist(values, bins=bins, color="seagreen", edgecolor="white")
    else:
        ax.text(
            0.5,
            0.5,
            "No usable objects",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )

    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.grid(axis="y", alpha=0.25)
    annotate_used_objects(ax, len(values), n_total)


def plot_highlighted_fraction_by_band(delta_metrics, output_dir):
    fig, axes = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(12, 8))
    axes = axes.flatten()
    n_total_objects = len(delta_metrics)
    n_used_by_band = {}

    for i, band in enumerate(BANDS):
        values = values_from_column(delta_metrics, f"highlighted_fraction_{band}")
        n_used_by_band[band] = len(values)
        plot_fraction_histogram(ax=axes[i], values=values, title=f"{band}-band", n_total=n_total_objects)
        axes[i].set_ylabel("Number of objects")
        if i >= 4:
            axes[i].set_xlabel("Highlighted / total LSST points")

    fig.suptitle("Highlighted-Point Fraction by Band")
    output_path = os.path.join(output_dir, "highlighted_fraction_by_band.png")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save_figure(fig, output_path)
    plt.close(fig)

    return {
        "highlighted_fraction_by_band_plot": output_path,
        "n_used_highlighted_fraction_by_band": n_used_by_band,
    }


def plot_global_highlighted_fraction(delta_metrics, output_dir):
    values = values_from_column(delta_metrics, "highlighted_fraction")

    fig, ax = plt.subplots(figsize=(8, 5))
    plot_fraction_histogram(
        ax=ax,
        values=values,
        title="Global Highlighted-Point Fraction",
        n_total=len(delta_metrics),
    )
    ax.set_xlabel("Highlighted / total LSST points")
    ax.set_ylabel("Number of objects")

    output_path = os.path.join(output_dir, "highlighted_fraction_global.png")
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close(fig)

    return {
        "highlighted_fraction_global_plot": output_path,
        "n_used_highlighted_fraction_global": len(values),
    }


def plot_count_histogram(ax, values, title):
    if values:
        ax.hist(values, bins=integer_bins(values), color="slateblue", edgecolor="white")
    else:
        ax.text(
            0.5,
            0.5,
            "No objects",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )

    ax.set_title(title)
    ax.set_xlim(left=-0.5)
    ax.grid(axis="y", alpha=0.25)


def plot_lsst_points_by_band(delta_metrics, output_dir):
    fig, axes = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(12, 8))
    axes = axes.flatten()
    n_used_by_band = {}

    for i, band in enumerate(BANDS):
        values = values_from_column(delta_metrics, f"n_total_points_{band}")
        n_used_by_band[band] = len(values)
        plot_count_histogram(axes[i], values, f"{band}-band")
        axes[i].set_ylabel("Number of objects")
        if i >= 4:
            axes[i].set_xlabel("Number of LSST points")

    fig.suptitle("Number of LSST Points by Band")
    output_path = os.path.join(output_dir, "lsst_points_by_band.png")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save_figure(fig, output_path)
    plt.close(fig)

    return {
        "lsst_points_by_band_plot": output_path,
        "n_used_lsst_points_by_band": n_used_by_band,
    }


def plot_global_lsst_points(delta_metrics, output_dir):
    values = values_from_column(delta_metrics, "n_total_points")
    positive_values = [value for value in values if value > 0]
    n_zero_or_negative = len(values) - len(positive_values)

    fig, ax = plt.subplots(figsize=(8, 5))
    if positive_values:
        ax.hist(
            positive_values,
            bins=logarithmic_bins(positive_values),
            color="slateblue",
            edgecolor="white",
        )
        ax.set_xscale("log")
    else:
        ax.text(
            0.5,
            0.5,
            "No positive LSST point counts",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )

    if n_zero_or_negative:
        ax.text(
            0.98,
            0.95,
            f"Non-positive counts excluded: {n_zero_or_negative}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
        )

    ax.set_title("Global Number of LSST Points")
    ax.set_xlabel("Number of LSST points")
    ax.set_ylabel("Number of objects")
    ax.grid(axis="y", alpha=0.25)

    output_path = os.path.join(output_dir, "lsst_points_global.png")
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close(fig)

    return {
        "lsst_points_global_plot": output_path,
        "n_used_lsst_points_global": len(values),
    }


def main(delta_metrics_csv, output_dir, output_subdir="pop_plots"):
    output_dir = os.path.join(output_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)

    delta_metrics = load_delta_metrics(delta_metrics_csv)
    result = plot_delta_distribution(delta_metrics, output_dir)
    result.update(plot_highlighted_fraction_by_band(delta_metrics, output_dir))
    result.update(plot_global_highlighted_fraction(delta_metrics, output_dir))
    result.update(plot_lsst_points_by_band(delta_metrics, output_dir))
    result.update(plot_global_lsst_points(delta_metrics, output_dir))
    result["pop_plots_dir"] = output_dir

    print(f"Saved delta distribution to {result['delta_distribution_plot']}")
    print(f"Saved highlighted-point fraction by band to {result['highlighted_fraction_by_band_plot']}")
    print(f"Saved global highlighted-point fraction to {result['highlighted_fraction_global_plot']}")
    print(f"Saved LSST point counts by band to {result['lsst_points_by_band_plot']}")
    print(f"Saved global LSST point counts to {result['lsst_points_global_plot']}")
    print(
        "Used objects for relative time range of interest: "
        f"{result['n_used_deltas']}/"
        f"{result['n_used_deltas'] + result['n_single_point_deltas']}"
    )

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create population-level BOOM run plots.")
    parser.add_argument("--delta_metrics_csv", required=True, help="Path to delta_metrics.csv from plot_lc.py")
    parser.add_argument("--output_dir", required=True, help="Run directory where pop_plots is created")
    parser.add_argument("--output_subdir", default="pop_plots", help="Population plot subfolder name")
    args = parser.parse_args()

    main(
        args.delta_metrics_csv,
        args.output_dir,
        output_subdir=args.output_subdir,
    )
