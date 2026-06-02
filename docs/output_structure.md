# Output Structure

The pipeline writes one self-contained run directory for each execution.

```text
boom_filter_outputs/
  run_<creation-time>_range_<start>_to_<end>/
    cands.json
    objs.json
    cands_sum.csv
    objs_sum.csv
    plots/
      lc_<n>_object_<objectId>.png
    pipeline_manifest.json
```

## Files

`cands.json`: full candidate records returned by the BOOM candidate query.

`objs.json`: full object records returned by the BOOM object query.

`cands_sum.csv`: compact table of candidate IDs, object IDs, detection time, magnitude, and SNR.

`objs_sum.csv`: compact table of object-level photometry counts, time baselines, coordinates, magnitudes, SNRs, and per-band counts.

`plots/`: light-curve PNGs. Each plot shows LSST detections, highlighted points that passed the filter, LSST forced photometry, and ZTF crossmatch detections when enabled and available.

`pipeline_manifest.json`: machine-readable record of the command options and artifact paths for the run.
