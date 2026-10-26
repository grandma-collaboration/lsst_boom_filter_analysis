# BOOM LSST Filter Analysis

This repository contains a small analysis pipeline for running BOOM LSST alert
filters, retrieving the matching candidates and objects, and producing
light-curve and population-level diagnostic plots.

The recommended entry point is `scripts/boom_pipeline.py`, configured through
`scripts/inputs.py`. In one run it:

1. runs a BOOM filter over a requested UTC time range, optionally skipping one
   excluded interval;
2. downloads candidate and object query products, and drops the objects with
   photometry in the excluded interval;
3. optionally counts every LSST alert and object in the range, to report the
   fraction passing the filter;
4. fetches BOOM's latest rising/fading rate in each band for each object;
5. writes compact CSV summaries;
6. computes per-object light-curve metrics;
7. optionally writes individual light-curve PNG and pickle files;
8. writes population-level plots and a markdown run summary.

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

Edit `scripts/inputs.py` (set at least `FILTER_FILE`, `START_TIME` and
`END_TIME`), then run the full pipeline from the repository root:

```bash
python scripts/boom_pipeline.py
```

`scripts/boom_pipeline.py` has no command-line interface: every parameter is
read from `scripts/inputs.py`, and any argument passed on the command line is
ignored.

With `INCLUDE_FORCED_PHOTOMETRY` and `INCLUDE_ZTF` set to `True`, as they are in
`scripts/inputs.py`, individual lightcurve plots include LSST detections, LSST
forced photometry, and ZTF cone-search crossmatch photometry when available. For
population-level diagnostic plots, only the LSST detections are included.

## Main Pipeline Parameters

`scripts/boom_pipeline.py` is the main command to use for most runs. All its
parameters are defined in `scripts/inputs.py` and must stay defined there. A
minimal configuration sets:

```python
FILTER_FILE = "configs/schecks.json"
START_TIME = "2026-03-01T12:00:00"
END_TIME = "2026-05-13T12:00:00"
```

The optional parameters can be left at their typical value.

| Parameter | Required | Typical value | Description |
| --- | --- | --- | --- |
| `FILTER_FILE` | yes | none | Path to a JSON file containing a MongoDB aggregation pipeline. The file must contain a JSON list. Relative paths are resolved from the current working directory. |
| `START_TIME` | yes | none | UTC ISO start time for the BOOM filter query, for example `"2026-03-01T12:00:00"`. |
| `END_TIME` | yes | none | UTC ISO end time for the BOOM filter query. It must be later than `START_TIME`. |
| `OUTPUT_DIR` | no | `"boom_filter_outputs"` | Base output directory. A timestamped run directory is created inside it. |
| `PLOTS_SUBDIR` | no | `"plots"` | Subdirectory inside the run directory where per-object plots and `delta_metrics.csv` are written. It must stay inside the run directory; absolute paths and `..` are rejected. |
| `EXCLUDE_START_TIME` | no | `None` | UTC ISO start time of an interval to exclude from the run (see [Time Range Exclusion](#time-range-exclusion)). Must be set together with `EXCLUDE_END_TIME`. |
| `EXCLUDE_END_TIME` | no | `None` | UTC ISO end time of the excluded interval. Must be set together with `EXCLUDE_START_TIME`. |
| `MAX_OBJECTS` | no | `None` | Limit the number of objects processed by the plotting and metric step. Useful for smoke tests. If set, it must be a positive integer. |
| `INCLUDE_ZTF` | no | `True` | Include ZTF cone-search crossmatch photometry in individual light-curve plots, when available. |
| `INCLUDE_FORCED_PHOTOMETRY` | no | `True` | Include LSST forced photometry in individual light-curve plots. |
| `SKIP_INDIVIDUAL_PLOTS` | no | `False` | Skip per-object light-curve PNG and pickle files. The pipeline still writes LSST-only `delta_metrics.csv`, population plots, query summary, and manifest. |
| `COMPUTE_FILTER_PASS_RATES` | no | `False` | Also count every LSST alert and object in the queried range, to report the fraction passing the filter (see [Filter Pass Rates](#filter-pass-rates)). Adds about 2 minutes. |

### Time Range Exclusion

Set `EXCLUDE_START_TIME` and `EXCLUDE_END_TIME` together in `scripts/inputs.py`
to exclude an interval inside the requested query range:

```python
EXCLUDE_START_TIME = "2026-03-05T12:00:00"
EXCLUDE_END_TIME = "2026-03-08T12:00:00"
```

The exclusion is applied in two steps:

1. The BOOM filter is only run outside the interval, so no alert inside it can
   pass the filter.
2. After the objects are retrieved, every object with any point inside
   `[EXCLUDE_START_TIME, EXCLUDE_END_TIME)` (LSST detection or forced
   photometry, in any band) is dropped from the run, along with its candidates.
   This includes forced-photometry measurements without a magnitude, which are
   never plotted. Objects with no point inside the interval, including objects
   with points on both sides of it, are kept with their full history.

Everything downstream (JSON and CSV files, light-curve metrics, plots and
summaries) only sees the objects that were kept. ZTF crossmatch photometry is
not affected by the exclusion. On individual light-curve plots, the interval is
shaded red where it overlaps the plotted time range.

The exclusion interval must have an end time after its start time. It cannot
remove the entire requested query interval. For now this is limited to only excluding one interval.

### Filter Pass Rates

With `COMPUTE_FILTER_PASS_RATES = True`, the pipeline also uses BOOM count
queries to count:

- every LSST alert detected in the queried range, minus the alerts inside the
  excluded interval;
- every LSST object with at least one detection in the queried range and no
  point (detection or forced photometry) inside the excluded interval, which is
  the same rule used to drop filter-passing objects.

It prints the number of filter-passing alerts and objects as fractions of these
totals, and saves them in `pipeline_manifest.json` under
`filter_result.filter_pass_rates`. The object count scans the whole
`LSST_alerts_aux` catalog, which takes about a minute.

To match the alert total, the filter-passing alert count only leaves out the
alerts inside the excluded interval (which the filter never sees), not the
alerts of dropped objects. It can therefore be larger than the number of
candidates in `cands.json`.

### Smoke Test Example

To run the full query but only process the first 10 objects during the plotting
and metric step, set in `scripts/inputs.py`:

```python
MAX_OBJECTS = 10
```

To quickly produce run-level metrics and population plots without individual
light-curve image files, set:

```python
SKIP_INDIVIDUAL_PLOTS = True
```

## Filter Configurations

An example filter pipeline is stored in `configs/`:

```text
configs/schecks.json
```

The .json file (`FILTER_FILE`) is passed to BOOM as the `pipeline` field of the filter request. The
pipeline runner queries BOOM in one-day batches and uses a request limit of
`10000` alerts per batch (`FILTER_QUERY_LIMIT` in
`scripts/boom_filter_per_run.py`). If a batch returns `10000` alerts, the script
prints a warning because the result may have been truncated by that limit.

## Outputs

Each full pipeline run creates one timestamped folder under `OUTPUT_DIR`:

```text
boom_filter_outputs/
  run_<creation-time>_range_<start>_to_<end>/
    cands.json
    objs.json
    cands_sum.csv
    objs_sum.csv
    band_rates.json
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
- `band_rates.json`: BOOM's rising and fading rates for each object and band, as
  of the object's most recent LSST detection in that band. They are shown in the
  legend of individual light-curve plots.
- `plots/delta_metrics.csv`: per-object metric table used to build the
  population plots.
- `pop_plots/`: run-level diagnostic plots.
- `query_summary.md`: markdown summary of retrieved candidates/objects and the
  number of objects with at least N filter-passing points.
- `pipeline_manifest.json`: machine-readable record of the parameters,
  generated artifact paths and, when enabled, filter pass rates.

When an excluded interval is set, all of these files only describe the objects
kept after the exclusion. See `docs/output_structure.md`.

## Standalone Scripts

The full pipeline calls the scripts below internally. Unlike
`boom_pipeline.py`, they take command-line arguments, and can be run separately
when you want to rerun only one step.

### `scripts/boom_filter_per_run.py`

Runs only the BOOM filter and writes candidate/object JSON, `band_rates.json`,
and CSV summaries.

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--filter_file` | yes | none | Path to the MongoDB pipeline JSON file. |
| `--start_time` | yes | none | ISO start time for the BOOM filter query. |
| `--end_time` | yes | none | ISO end time for the BOOM filter query. |
| `--output_dir` | no | `boom_filter_outputs` | Base output directory. A timestamped run folder is created inside it. |
| `--exclude_start_time` | no | `None` | Optional ISO start time of the excluded interval. It is not queried, and objects with any point in it are dropped. |
| `--exclude_end_time` | no | `None` | Optional ISO end time of the excluded interval. |
| `--compute_pass_rates` | no | `False` | Also count every LSST alert and object in the range and report the filter pass fractions. |

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
| `--exclude_start_time` | no | `None` | Optional UTC ISO start time of the run's excluded interval, only used to shade it on plots. Objects with points in it are dropped by `boom_filter_per_run.py`, not here. |
| `--exclude_end_time` | no | `None` | Optional UTC ISO end time of the excluded interval. |
| `--band_rates_json` | no | `None` | Optional path to `band_rates.json`, to show each band's rising/fading rate in the plot legends. |

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
- The full pipeline includes ZTF photometry according to `INCLUDE_ZTF` in
  `scripts/inputs.py`, while standalone `plot_lc.py` only includes it with
  `--include_ztf`.
- Standalone `plot_lc.py` does not read the excluded interval or the band rates
  from the run: pass `--exclude_start_time`/`--exclude_end_time` and
  `--band_rates_json` to get the shading and the rate legend entries.
