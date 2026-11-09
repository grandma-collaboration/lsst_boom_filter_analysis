"""
Input file for scripts/ghost/boom_filter_ghost.py.

Same layout as scripts/inputs.py: every parameter is set here, and the
script has no command-line options. Edit the values below, then, from the
root of the repository, run:

    python3 scripts/ghost/boom_filter_ghost.py

BOOM_USERNAME and BOOM_PASSWORD must be set in the shell, as for
scripts/boom_pipeline.py.

Extra packages beyond requirements.txt: pandas (always); scipy, astroquery
and astro_ghost (host association); astroquery (INSPECT). After installing
astro_ghost, run scripts/ghost/patch_ghost_dlr.py once: it fixes two bugs in
astro_ghost's DLR.py.
"""

from __future__ import annotations

# =========================================================================
# REQUIRED PARAMETERS
# =========================================================================

# Path to the JSON file containing the MongoDB aggregation pipeline that
# defines the BOOM filter. Must be a JSON list, same format as configs/*.json.
# Relative paths are resolved from the current working directory.
# A private filter should live outside this repository.
FILTER_FILE = "configs/schecks.json"

# UTC ISO start time of the BOOM query. The filter is run in 1-day batches.
START_TIME = "2026-02-26T12:00:00"

# UTC ISO end time of the BOOM query. Must be later than START_TIME.
END_TIME = "2026-03-08T12:00:00"

# =========================================================================
# OPTIONAL PARAMETERS
# =========================================================================

# Base output directory. The script writes into four subfolders of it:
# candidate_quicklooks/ (candidate CSVs), ned_flagged/ (candidates with an
# associated host and a usable distance modulus, plus the INSPECT HTML),
# candidate_jsons/ and object_jsons/ (only with SAVE_JSON). File names carry
# the run time and the queried range, so runs don't overwrite each other.
OUTPUT_DIR = "boom_ghost_outputs"

# If True, also download the objects and save the full candidate and object
# JSON.
SAVE_JSON = False

# If True, skip host association (GHOST, NED-LVS, Legacy DR10): faster, the
# host columns of the candidate CSVs stay empty, and nothing is written to
# ned_flagged/ (INSPECT is ignored).
SKIP_HOST_ASSOCIATION = False

# If True, re-query NED and SIMBAD around each flagged candidate and write an
# HTML page with Aladin views. Runs only with host association.
INSPECT = False

# Cone radius (arcsec) for the INSPECT NED/SIMBAD queries.
INSPECT_RADIUS_ARCSEC = 60.0

# Aladin field of view (arcmin) on the INSPECT page.
INSPECT_FOV_ARCMIN = 3.0
