# BOOM LSST Filter Analysis

This repository contains a small analysis pipeline for running BOOM LSST alert
filters, retrieving the matching candidates and objects, and producing
light-curve and population-level diagnostic plots.

The recommended entry point is `scripts/boom_pipeline.py`. In one command it:

1. runs a BOOM filter over a requested UTC time range;
2. downloads candidate and object query products;
3. writes compact CSV summaries;
4. computes per-object light-curve metrics;
5. optionally writes individual light-curve PNG and pickle files;
6. writes population-level plots and a markdown run summary.

BOOM credentials are required and are read from environment variables.

## Setup

Create and activate a Python environment, then install the dependencies:

```bash
conda create -n boom-queries python=3.11
conda activate boom-queries
pip install -r requirements.txt
```

Set your BOOM credentials before running any script that queries BOOM:

```bash
export BOOM_USERNAME="your_boom_username_here"
export BOOM_PASSWORD="your_boom_password_here"
```

You can also copy `.env.example` as a local template, but the scripts do not
load `.env` automatically. The variables must be present in the shell
environment.

## Quick Start

Run the full pipeline from the repository root:

```bash
python scripts/boom_pipeline.py \
  --output_dir boom_filter_outputs \
  --filter_file configs/schecks.json \
  --start_time 2026-03-01T12:00:00 \
  --end_time 2026-05-13T12:00:00
```

By default, for individual lightcurve plots, the full pipeline includes LSST detections, LSST forced photometry,
and ZTF cone-search crossmatch photometry when available. For population-level diagnostic plots, only the LSST detections are included.

## Main Pipeline Arguments

`scripts/boom_pipeline.py` is the main command to use for most runs.

```bash
python scripts/boom_pipeline.py \
  --filter_file configs/schecks.json \
  --start_time 2026-03-01T12:00:00 \
  --end_time 2026-05-13T12:00:00
```

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--filter_file` | yes | none | Path to a JSON file containing a MongoDB aggregation pipeline. The file must contain a JSON list. |
| `--start_time` | yes | none | UTC ISO start time for the BOOM filter query, for example `2026-03-01T12:00:00`. |
| `--end_time` | yes | none | UTC ISO end time for the BOOM filter query. It must be later than `--start_time`. |
| `--output_dir` | no | `boom_filter_outputs` | Base output directory. A timestamped run directory is created inside it. |
| `--exclude_start_time` | no | `None` | UTC ISO start time of an interval to exclude from the BOOM query. Must be used with `--exclude_end_time`. |
| `--exclude_end_time` | no | `None` | UTC ISO end time of an interval to exclude from the BOOM query. Must be used with `--exclude_start_time`. |
| `--plots_subdir` | no | `plots` | Subdirectory inside the run directory where per-object plots and `delta_metrics.csv` are written. It must stay inside the run directory; absolute paths and `..` are rejected. |
| `--no_ztf` | no | `False` | Disable ZTF cone-search crossmatch photometry in individual light-curve plots. Without this flag, ZTF photometry is included when available. |
| `--no_forced_photometry` | no | `False` | Disable LSST forced photometry in individual light-curve plots. Without this flag, forced photometry is included. |
| `--max_objects` | no | `None` | Limit the number of objects processed by the plotting and metric step. Useful for smoke tests. If provided, it must be a positive integer. |
| `--skip_individual_plots` | no | `False` | Skip per-object light-curve PNG and pickle files. The pipeline still writes LSST-only `delta_metrics.csv`, population plots, query summary, and manifest. |

### Time Range Exclusion

Use `--exclude_start_time` and `--exclude_end_time` together to skip an interval
inside the requested query range:

```bash
python scripts/boom_pipeline.py \
  --filter_file configs/schecks.json \
  --start_time 2026-03-01T12:00:00 \
  --end_time 2026-05-13T12:00:00 \
  --exclude_start_time 2026-03-05T12:00:00 \
  --exclude_end_time 2026-03-08T12:00:00
```

The exclusion interval must have an end time after its start time. It cannot
remove the entire requested query interval. For now this is limited to only excluding one interval.

### Smoke Test Example

This command runs the full query but only processes the first 10 objects during
the plotting and metric step:

```bash
python scripts/boom_pipeline.py \
  --filter_file configs/schecks.json \
  --start_time 2026-03-01T12:00:00 \
  --end_time 2026-05-13T12:00:00 \
  --max_objects 10
```

To quickly produce run-level metrics and population plots without individual
light-curve image files:

```bash
python scripts/boom_pipeline.py \
  --filter_file configs/schecks.json \
  --start_time 2026-03-01T12:00:00 \
  --end_time 2026-05-13T12:00:00 \
  --skip_individual_plots
```

## Filter Configurations

An example filter pipeline is stored in `configs/`:

```text
configs/schecks.json
```

The .json file (--filter_file) is passed to BOOM as the `pipeline` field of the filter request. The
pipeline runner queries BOOM in one-day batches and uses a request limit of
`10000` alerts per batch. If a batch returns `10000` alerts, the script prints a
warning because the result has been truncated by that limit.

## Outputs

Each full pipeline run creates one timestamped folder under `--output_dir`:

```text
boom_filter_outputs/
  run_<creation-time>_range_<start>_to_<end>/
    cands.json
    objs.json
    cands_sum.csv
    objs_sum.csv
    plots/
      lc_<n>_object_<objectId>.png
      lc_<n>_object_<objectId>.pickle
      delta_metrics.csv
    pop_plots/
      delta_distribution.png
      delta_distribution.pickle
      highlighted_fraction_by_band.png
      highlighted_fraction_by_band.pickle
      highlighted_fraction_global.png
      highlighted_fraction_global.pickle
      lsst_points_by_band.png
      lsst_points_by_band.pickle
      lsst_points_global.png
      lsst_points_global.pickle
    query_summary.md
    pipeline_manifest.json
```

Important files:

- `cands.json`: full candidate records returned by BOOM.
- `objs.json`: full object records returned by BOOM.
- `cands_sum.csv`: compact candidate table with candidate ID, object ID, JD,
  magnitude, and SNR.
- `objs_sum.csv`: compact object table with photometry counts, time baselines,
  coordinates, magnitudes, SNRs, and per-band counts.
- `plots/delta_metrics.csv`: per-object metric table used to build the
  population plots.
- `pop_plots/`: run-level diagnostic plots.
- `query_summary.md`: markdown summary of retrieved candidates/objects and the
  number of objects with at least N filter-passing points.
- `pipeline_manifest.json`: machine-readable record of command options and
  generated artifact paths.

See `docs/output_structure.md`.

## Standalone Scripts

The full pipeline calls the scripts below internally. They can also be run
separately when you want to rerun only one step.

### `scripts/boom_filter_per_run.py`

Runs only the BOOM filter and writes candidate/object JSON plus CSV summaries.

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--filter_file` | yes | none | Path to the MongoDB pipeline JSON file. |
| `--start_time` | yes | none | ISO start time for the BOOM filter query. |
| `--end_time` | yes | none | ISO end time for the BOOM filter query. |
| `--output_dir` | no | `boom_filter_outputs` | Base output directory. A timestamped run folder is created inside it. |
| `--exclude_start_time` | no | `None` | Optional ISO start time of an excluded query interval. |
| `--exclude_end_time` | no | `None` | Optional ISO end time of an excluded query interval. |

### `scripts/plot_lc.py`

Plots light curves and writes `delta_metrics.csv` from existing `objs.json` and
`cands.json` files.

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--objects_json` | yes | none | Path to an objects JSON file, normally `objs.json`. |
| `--alerts_json` | yes | none | Path to a candidate/alert JSON file, normally `cands.json`. |
| `--output_dir` | no | `./lc_plots/` | Directory where the plot subdirectory is created. |
| `--include_ztf` | no | `False` | Include ZTF cone-search crossmatch photometry when available. |
| `--no_fp` | no | `False` | Disable LSST forced photometry. By default forced photometry is included. |
| `--output_subdir` | no | `None` | Optional plot subfolder name. If omitted, it is derived from the objects JSON filename. |
| `--max_objects` | no | `None` | Optional maximum number of objects to process. |
| `--skip_individual_plots` | no | `False` | Only write LSST-only delta metrics; do not write individual light-curve plots. |

### `scripts/population_plots.py`

Creates population-level plots from an existing `delta_metrics.csv`.

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--delta_metrics_csv` | yes | none | Path to `delta_metrics.csv` from `plot_lc.py` or the full pipeline. |
| `--output_dir` | yes | none | Run directory where the population plot subdirectory is created. |
| `--output_subdir` | no | `pop_plots` | Population plot subfolder name. |

### `scripts/run_summary.py`

Writes a markdown summary for one run.

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--objects_summary_csv` | yes | none | Path to `objs_sum.csv`. |
| `--output_dir` | yes | none | Run directory where the markdown summary is written. |
| `--n_candidates` | yes | none | Total number of candidates retrieved. |
| `--n_objects` | yes | none | Total number of objects retrieved. |
| `--output_filename` | no | `query_summary.md` | Name of the summary markdown file. |
| `--delta_metrics_csv` | no | `None` | Optional path to `plots/delta_metrics.csv`; enables the first-relative-time-range-bin section. |

## Notes

- Do not commit BOOM credentials or generated run outputs.
- `--no_ztf` on the full pipeline and `--include_ztf` on standalone
  `plot_lc.py` differ because the full pipeline opts into ZTF by default, while
  standalone plotting does not.
