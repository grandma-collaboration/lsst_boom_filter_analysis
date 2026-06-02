# BOOM Queries And Plots

Small pipeline for running BOOM LSST alert filters, saving candidate/object query products, and producing light-curve plots with LSST detections, LSST forced photometry, and ZTF crossmatch photometry. BOOM credentials are needed, setup as environment variables.

## Setup

Create and activate a conda environment, then install the required dependencies:

```bash
conda create -n boom-queries
conda activate boom-queries
pip install -r requirements.txt
export BOOM_USERNAME="your_boom_username_here"
export BOOM_PASSWORD="your_boom_password_here"
```

## Run

```bash
python scripts/boom_pipeline.py \
  --output_dir boom_filter_outputs \
  --filter_file configs/.json \
  --start_time 2026-03-01T12:00:00 \
  --end_time 2026-05-13T12:00:00
```

By default, the plots include LSST detections, LSST forced photometry, and ZTF crossmatch photometry when available.

Useful options:

```bash
--no_ztf
--no_forced_photometry
--max_objects 10
```

## Outputs

Each run creates one timestamped folder under `--output_dir`:

```text
run_.../
  cands.json
  objs.json
  cands_sum.csv
  objs_sum.csv
  plots/
  pipeline_manifest.json
```

See [docs/output_structure.md](docs/output_structure.md) for details.

## Notes

Do not commit credentials or generated run outputs. Use `.env.example` as a template for local environment setup.
