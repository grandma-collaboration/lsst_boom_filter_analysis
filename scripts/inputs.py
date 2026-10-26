"""
Input file for scripts/boom_pipeline.py.

Every parameter the pipeline accepts is defined here, so this file is the
single place to edit before a run. boom_pipeline.py has no command-line
interface of its own and reads every value straight from this module.

Edit the values below, then, from the root of the repository, run:

    python3 scripts/boom_pipeline.py
"""

from __future__ import annotations

# =========================================================================
# REQUIRED PARAMETERS
# =========================================================================
# The pipeline cannot run without these being set to something valid.

# Path to the JSON file containing the MongoDB aggregation pipeline that
# defines the BOOM filter to test. Must be a JSON list. Relative paths are
# resolved from the current working directory, not from this file.
FILTER_FILE = "configs/schecks.json"

# UTC ISO start time of the BOOM query, e.g. "2026-03-01T12:00:00".
START_TIME = "2025-09-01T12:00:00"
# START_TIME = "2026-07-06T12:00:00"

# UTC ISO end time of the BOOM query. Must be strictly later than START_TIME.
END_TIME = "2026-07-15T12:00:00"
# END_TIME = "2026-07-11T12:00:00"

# =========================================================================
# OPTIONAL PARAMETERS
# =========================================================================
# These all have a meaningful "off" value (None or False) and don't have to
# be changed.

# Base output directory. A new timestamped run folder is created inside it
# for every run, so this can be reused across runs without overwriting past
# results.
OUTPUT_DIR = "boom_filter_outputs"

# Subdirectory inside the run folder where individual light-curve plots and
# delta_metrics.csv are written. Must stay inside the run directory: absolute
# paths and ".." are rejected.
PLOTS_SUBDIR = "plots"

# Optional UTC ISO start/end time of one interval to exclude from the run.
# Both must be set together, or both left as None to disable exclusion.
# In this version of the pipeline, an object with any point -- any band,
# detection or forced photometry -- inside [EXCLUDE_START_TIME,
# EXCLUDE_END_TIME) is dropped from the run entirely, as if it had never been
# retrieved. Objects that don't touch this range at all are kept with their
# full, untouched history. The excluded range is still shaded red on plots
# for reference, even though no retrieved object ever has data in it.
EXCLUDE_START_TIME = "2026-02-18T12:00:00"
EXCLUDE_END_TIME = "2026-02-27T12:00:00"
# EXCLUDE_START_TIME = None
# EXCLUDE_END_TIME = None

# Optional cap on the number of objects processed by the plotting and metric
# step, useful for smoke tests. Set to None to process every object returned
# by the filter. If set, must be a positive integer.
# MAX_OBJECTS = 25
MAX_OBJECTS = None

# Whether to include ZTF cone-search crossmatch photometry in individual
# light-curve plots, when a match is available. This is the only place the
# exclusion range above does not apply: ZTF points are never dropped by it.
INCLUDE_ZTF = True

# Whether to include LSST forced photometry in individual light-curve plots.
INCLUDE_FORCED_PHOTOMETRY = True

# If True, skip writing individual light-curve PNG/pickle files entirely.
# The pipeline still writes LSST-only delta_metrics.csv, population plots,
# the query summary, and the manifest. Useful for large runs where only the
# population-level outputs are needed.
SKIP_INDIVIDUAL_PLOTS = False

# If True, also count every LSST alert and every LSST object in the queried
# range, and report the fraction of them that pass the filter (printed after
# the filter queries, and saved as "filter_pass_rates" in the manifest). The
# excluded range above applies to these totals the same way as to the
# filter-passing objects: alerts inside it are not counted, and objects with
# any point (detection or forced photometry) inside it are not counted.
COMPUTE_FILTER_PASS_RATES = True
